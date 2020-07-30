/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#undef SHOW_WARNINGS							// Display compilation warnings
#undef USE_DOUBLE								// If enabled, double-precesion is used
#undef FAST_COMPILE								// If enabled, only a single version of the reconstruction code is compiled
#undef ARRAY_DEBUG								// If enabled, array access is tested for validity
#define DATA_DEGREE 0							// The order of the B-Spline used to splat in data for color interpolation
												// This can be changed to zero if more interpolatory performance is desired.
#define WEIGHT_DEGREE 2							// The order of the B-Spline used to splat in the weights for density estimation
#define NORMAL_DEGREE 2							// The order of the B-Spline used to splat in the normals for constructing the Laplacian constraints
#define DEFAULT_FEM_DEGREE 1					// The default finite-element degree
#define DEFAULT_FEM_BOUNDARY BOUNDARY_NEUMANN	// The default finite-element boundary type
#define DIMENSION 3								// The dimension of the system

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "PPolynomial.h"
#include "FEMTree.h"
#include "Ply.h"
#include "PointStreamData.h"

#include "../AdaptiveSolvers.h"
#include "../my_helpers.h"

MessageWriter messageWriter;

const float DefaultPointWeightMultiplier = 2.f;

cmdLineParameter< char* >
	In( "in" ) ,
	Out( "out" ) ,
	TempDir( "tempDir" ) ,
	VoxelGrid( "voxel" ) ,
	Tree( "tree" ) ,
	Transform( "xForm" );

cmdLineReadable
	Performance( "performance" ) ,
	ShowResidual( "showResidual" ) ,
	NoComments( "noComments" ) ,
	PolygonMesh( "polygonMesh" ) ,
	NonManifold( "nonManifold" ) ,
	ASCII( "ascii" ) ,
	Density( "density" ) ,
	LinearFit( "linearFit" ) ,
	PrimalVoxel( "primalVoxel" ) ,
	ExactInterpolation( "exact" ) ,
	Normals( "normals" ) ,
	Colors("colors"),
	InCore("inCore"),
	LevelSetMap("levelSetMap"),
	SkipContouring("skipContouring"),
	Verbose("verbose");

cmdLineParameter< int >
#ifndef FAST_COMPILE
	Degree( "degree" , DEFAULT_FEM_DEGREE ) ,
#endif // !FAST_COMPILE
	Depth( "depth" , 8 ) ,
	KernelDepth( "kernelDepth", 6 ) ,
	FullDepth( "fullDepth" , 5 ) ,
	CGDepth("cgDepth", 0 ) ,
	BaseDepth("baseDepth", 0),
	Iters("iters", 8),
	BaseVCycles( "baseVCycles" , 1 ) ,
#ifndef FAST_COMPILE
	BType( "bType" , DEFAULT_FEM_BOUNDARY+1 ) ,
#endif // !FAST_COMPILE
	MaxMemoryGB( "maxMemory" , 0 ) ,
	Threads( "threads" , omp_get_num_procs() );

cmdLineParameter< float >
	DataX( "data" , 32.f ) ,
	SamplesPerNode( "samplesPerNode" , 1.5f ) ,
	Scale( "scale" , 1.1f ) ,
	Width( "width" , 0.f ) ,
	Confidence( "confidence" , 0.f ) ,
	ConfidenceBias( "confidenceBias" , 0.f ) ,
	CGSolverAccuracy( "cgAccuracy" , 1e-3f ) ,
	PointWeight( "pointWeight" , DefaultPointWeightMultiplier);

cmdLineReadable* params[] =
{
#ifndef FAST_COMPILE
	&Degree , &BType ,
#endif // !FAST_COMPILE
	&In , &Depth , &Out , &Transform ,
	&Width ,
	&Scale , &Verbose , &CGSolverAccuracy , &NoComments ,
	&KernelDepth , &SamplesPerNode , &Confidence , &NonManifold , &PolygonMesh , &ASCII , &ShowResidual ,
	&ConfidenceBias ,
	&BaseDepth , &BaseVCycles ,
	&CGDepth ,
	&PointWeight ,
	&VoxelGrid , &Threads ,
	&Tree ,
	&Density ,
	&FullDepth ,
	&Iters ,
	&DataX ,
	&Colors ,
	&Normals ,
	&LinearFit ,
	&PrimalVoxel ,
	&TempDir ,
	&ExactInterpolation ,
	&Performance ,
	&MaxMemoryGB ,
	&InCore ,
	&LevelSetMap ,
	&SkipContouring,
	NULL
};

template <typename __Coord>
struct MyInOutStream
{
	const std::vector<__Coord>& pos_pts;
	const std::vector<__Coord>& nrl_pts;
	const std::vector<__Coord>& pos_aux_pts;
	const std::vector<__Coord>& nrl_aux_pts;

	__ProcBuffers<__Coord> proc_buffers;

	MyInOutStream(const std::vector<__Coord>& _pos_pts, const std::vector<__Coord>& _nrl_pts,
		const std::vector<__Coord>& _pos_aux_pts, const std::vector<__Coord>& _nrl_aux_pts) : pos_pts(_pos_pts), nrl_pts(_nrl_pts), pos_aux_pts(_pos_aux_pts), nrl_aux_pts(_nrl_aux_pts) {}
};

void ShowUsage(char* ex)
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input points>] = %s\n" , In.name, In.value );
	printf( "\t[--%s <ouput triangle mesh>] = %s\n" , Out.name, Out.value );
	printf( "\t[--%s <ouput voxel grid>] = %s\n" , VoxelGrid.name, VoxelGrid.value);
	printf( "\t[--%s <ouput fem tree>] = %s\n" , Tree.name, Tree.value);
#ifndef FAST_COMPILE
	printf( "\t[--%s <b-spline degree>] = %d]\n" , Degree.name , Degree.value );
	printf( "\t[--%s <boundary type>] = %d]\n" , BType.name , BType.value );
	for( int i=0 ; i<BOUNDARY_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i+1 , BoundaryNames[i] );
#endif // !FAST_COMPILE
	printf( "\t========== DEPTH ==========\n");
	printf( "\t[--%s <maximum reconstruction depth>] = %d\n" , Depth.name , Depth.value );
	printf( "\t[--%s <full (computation) depth (this <= depth)>] = %d\n" , FullDepth.name , FullDepth.value );
	printf( "\t[--%s <kernel depth (this <= depth)> =] %d\n" , KernelDepth.name , KernelDepth.value );
	printf( "\t[--%s <coarse MG solver depth (this <== fullDepth)>] = %d\n" , BaseDepth.name , BaseDepth.value );
	printf(" \t[--%s <CG depth (~this use CG solver)>] = %d\n", CGDepth.name, CGDepth.value);
	printf( "\t==========================\n");
	printf( "\t[--%s <voxel width>] = %f\n" , Width.name, Width.value);
	printf( "\t[--%s <coarse MG solver v-cycles>] = %d\n" , BaseVCycles.name , BaseVCycles.value );
	printf( "\t[--%s <scale factor>] = %f\n" , Scale.name , Scale.value );
	printf( "\t[--%s <minimum number of samples per node>] = %f\n" , SamplesPerNode.name, SamplesPerNode.value );
	printf( "\t[--%s <interpolation weight> =] %f * <b-spline degree>\n" , PointWeight.name , PointWeight.value);
	printf( "\t[--%s <iterations>] = %d\n" , Iters.name , Iters.value );
	printf( "\t[--%s] = %s\n" , ExactInterpolation.name , ExactInterpolation.set? "true" : "false");
	printf( "\t[--%s <pull factor>] = %f\n" , DataX.name , DataX.value );
	printf( "\t[--%s] = %s\n" , Colors.name, Colors.set ? "true" : "false");
	printf( "\t[--%s] = %s\n" , Normals.name, Normals.set ? "true" : "false");
#ifdef _OPENMP
	printf( "\t[--%s <num threads> = %d]\n" , Threads.name , Threads.value );
#endif // _OPENMP
	printf( "\t[--%s <normal confidence exponent> = %f]\n" , Confidence.name , Confidence.value );
	printf( "\t[--%s <normal confidence bias exponent> = %f]\n" , ConfidenceBias.name , ConfidenceBias.value );
	printf( "\t[--%s] = %s\n" , NonManifold.name , NonManifold.set ? "true" : "false");
	printf( "\t[--%s] = %s\n" , PolygonMesh.name , PolygonMesh.set ? "true" : "false");
	printf( "\t[--%s <cg solver accuracy> = %f]\n" , CGSolverAccuracy.name , CGSolverAccuracy.value );
	printf( "\t[--%s <maximum memory (in GB)> = %d]\n" , MaxMemoryGB.name , MaxMemoryGB.value );
	printf( "\t[--%s] = %s\n" , Performance.name , Performance.set ? "true" : "false");
	printf( "\t[--%s] = %s\n" , Density.name , Density.set ? "true" : "false");
	printf( "\t[--%s] = %s\n" , LinearFit.name , LinearFit.set ? "true" : "false");
	printf( "\t[--%s] = %s\n" , PrimalVoxel.name , PrimalVoxel.set ? "true" : "false");
	printf( "\t[--%s] = %s\n" , ASCII.name , ASCII.set ? "true" : "false");
	printf( "\t[--%s] = %s\n" , NoComments.name , NoComments.set ? "true" : "false");
	printf( "\t[--%s] = %s\n" , TempDir.name , TempDir.set ? "true" : "false");
	printf( "\t[--%s] = %s\n" , InCore.name , InCore.set ? "true" : "false");
	printf( "\t[--%s] = %s\n" , Verbose.name , Verbose.set ? "true" : "false");
}

double Weight( double v , double start , double end )
{
	v = ( v - start ) / ( end - start );
	if     ( v<0 ) return 1.;
	else if( v>1 ) return 0.;
	else
	{
		// P(x) = a x^3 + b x^2 + c x + d
		//		P (0) = 1 , P (1) = 0 , P'(0) = 0 , P'(1) = 0
		// =>	d = 1 , a + b + c + d = 0 , c = 0 , 3a + 2b + c = 0
		// =>	c = 0 , d = 1 , a + b = -1 , 3a + 2b = 0
		// =>	a = 2 , b = -3 , c = 0 , d = 1
		// =>	P(x) = 2 x^3 - 3 x^2 + 1
		return 2. * v * v * v - 3. * v * v + 1.;
	}
}

template< unsigned int Dim , class Real >
struct FEMTreeProfiler
{
	FEMTree< Dim , Real >& tree;
	double t;

	FEMTreeProfiler( FEMTree< Dim , Real >& t ) : tree(t) { ; }
	void start( void ){ t = Time() , FEMTree< Dim , Real >::ResetLocalMemoryUsage(); }
	void print( const char* header ) const
	{
		FEMTree< Dim , Real >::MemoryUsage();
		if( header ) printf( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , (double)MemoryInfo::PeakMemoryUsageMB() );
		else         printf(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , (double)MemoryInfo::PeakMemoryUsageMB() );
	}
	void dumpOutput( const char* header ) const
	{
		FEMTree< Dim , Real >::MemoryUsage();
		if( header ) messageWriter( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , (double)MemoryInfo::PeakMemoryUsageMB() );
		else         messageWriter(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , (double)MemoryInfo::PeakMemoryUsageMB() );
	}
	void dumpOutput2( std::vector< std::string >& comments , const char* header ) const
	{
		FEMTree< Dim , Real >::MemoryUsage();
		if( header ) messageWriter( comments , "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , (double)MemoryInfo::PeakMemoryUsageMB() );
		else         messageWriter( comments ,    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , (double)MemoryInfo::PeakMemoryUsageMB() );
	}
};

template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetBoundingBoxXForm( Point< Real , Dim > min , Point< Real , Dim > max , Real scaleFactor )
{
	Point< Real , Dim > center = ( max + min ) / 2;
	Real scale = max[0] - min[0];
	for( int d=1 ; d<Dim ; d++ ) scale = std::max< Real >( scale , max[d]-min[d] );
	scale *= scaleFactor;
	for( int i=0 ; i<Dim ; i++ ) center[i] -= scale/2;
	XForm< Real , Dim+1 > tXForm = XForm< Real , Dim+1 >::Identity() , sXForm = XForm< Real , Dim+1 >::Identity();
	for( int i=0 ; i<Dim ; i++ ) sXForm(i,i) = (Real)(1./scale ) , tXForm(Dim,i) = -center[i];
	return sXForm * tXForm;
}
template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetBoundingBoxXForm( Point< Real , Dim > min , Point< Real , Dim > max , Real width , Real scaleFactor , int& depth )
{
	// Get the target resolution (along the largest dimension)
	Real resolution = ( max[0]-min[0] ) / width;
	for( int d=1 ; d<Dim ; d++ ) resolution = std::max< Real >( resolution , ( max[d]-min[d] ) / width );
	resolution *= scaleFactor;
	depth = 0;
	while( (1<<depth)<resolution ) depth++;

	Point< Real , Dim > center = ( max + min ) / 2;
	Real scale = (1<<depth) * width;

	for( int i=0 ; i<Dim ; i++ ) center[i] -= scale/2;
	XForm< Real , Dim+1 > tXForm = XForm< Real , Dim+1 >::Identity() , sXForm = XForm< Real , Dim+1 >::Identity();
	for( int i=0 ; i<Dim ; i++ ) sXForm(i,i) = (Real)(1./scale ) , tXForm(Dim,i) = -center[i];
	return sXForm * tXForm;
}

template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetPointXForm( InputPointStream< Real , Dim >& stream , Real width , Real scaleFactor , int& depth )
{
	Point< Real , Dim > min , max;
	stream.boundingBox( min , max );
	return GetBoundingBoxXForm( min , max , width , scaleFactor , depth );
}
template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetPointXForm( InputPointStream< Real , Dim >& stream , Real scaleFactor )
{
	Point< Real , Dim > min , max;
	stream.boundingBox( min , max );
	return GetBoundingBoxXForm( min , max , scaleFactor );
}

template< unsigned int Dim , typename Real >
struct ConstraintDual
{
	Real target , weight;
	ConstraintDual(Real t, Real w) : target(t), weight(w) { }
	CumulativeDerivativeValues< Real, Dim, 0 > operator()(const Point< Real, Dim >& p) const {
		// p 가 여기선 안 쓰인다... -_-??
		// p 에 대해, 나름의 weight 를 추가할 수 있다.
		// 이건 나중에...
		// dojo to do
		auto ret = CumulativeDerivativeValues< Real, Dim, 0 >(target*weight);  // Point< Real , CumulativeDerivatives< Dim , D >::Size >;
		return ret;
	};
};
template< unsigned int Dim , typename Real >
struct SystemDual
{
	Real weight;
	SystemDual( Real w ) : weight(w){ }
	CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim >& p , const CumulativeDerivativeValues< Real , Dim , 0 >& dValues ) const { return dValues * weight; };
	CumulativeDerivativeValues< double , Dim , 0 > operator()( const Point< Real , Dim >& p , const CumulativeDerivativeValues< double , Dim , 0 >& dValues ) const { return dValues * weight; };
};
template< unsigned int Dim >
struct SystemDual< Dim , double >
{
	typedef double Real;
	Real weight;
	SystemDual( Real w ) : weight(w){ }
	CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim >& p , const CumulativeDerivativeValues< Real , Dim , 0 >& dValues ) const { return dValues * weight; };
};

template< unsigned int Dim, typename Vertex , typename Real , unsigned int ... FEMSigs , typename ... SampleData , typename __Coord >
void ExtractMesh( MyInOutStream<__Coord>& my_io,  UIntPack< FEMSigs ... > , std::tuple< SampleData ... > , FEMTree< sizeof ... ( FEMSigs ) , Real >& tree , const DenseNodeData< Real , UIntPack< FEMSigs ... > >& solution , Real isoValue , const std::vector< typename FEMTree< sizeof ... ( FEMSigs ) , Real >::PointSample >* samples , std::vector< MultiPointStreamData< Real , PointStreamNormal< Real , Dim > , MultiPointStreamData< Real , SampleData ... > > >* sampleData , const typename FEMTree< sizeof ... ( FEMSigs ) , Real >::template DensityEstimator< WEIGHT_DEGREE >* density ,
	std::function< void ( Vertex& , Point< Real , Dim > , Real , MultiPointStreamData< Real , PointStreamNormal< Real , Dim > , MultiPointStreamData< Real , SampleData ... > > ) > SetVertex , std::vector< std::string > &comments , XForm< Real , sizeof...(FEMSigs)+1 > iXForm ,
	std::function< void ( float*, float*, Vertex& ) > SetOutBuf)
{
	//static const int Dim = sizeof ... ( FEMSigs );
	typedef UIntPack< FEMSigs ... > Sigs;
	typedef PointStreamNormal< Real , Dim > NormalPointSampleData;
	typedef MultiPointStreamData< Real , SampleData ... > AdditionalPointSampleData;
	typedef MultiPointStreamData< Real , NormalPointSampleData , AdditionalPointSampleData > TotalPointSampleData;
	static const unsigned int DataSig = FEMDegreeAndBType< DATA_DEGREE , BOUNDARY_FREE >::Signature;
	typedef typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;

	FEMTreeProfiler< Dim , Real > profiler( tree );

	char tempHeader[1024];
	{
		char tempPath[1024];
		tempPath[0] = 0;
		if( TempDir.set ) strcpy( tempPath , TempDir.value );
		else SetTempDirectory( tempPath , sizeof(tempPath) );
		if( strlen(tempPath)==0 ) sprintf( tempPath , ".%c" , FileSeparator );
		if( tempPath[ strlen( tempPath )-1 ]==FileSeparator ) sprintf( tempHeader , "%sPR_" , tempPath );
		else                                                  sprintf( tempHeader , "%s%cPR_" , tempPath , FileSeparator );
	}
	CoredMeshData< Vertex > *mesh;
	if( InCore.set ) mesh = new CoredVectorMeshData< Vertex >();
	else             mesh = new CoredFileMeshData< Vertex >( tempHeader );

	profiler.start();
	typename IsoSurfaceExtractor< Dim , Real , Vertex >::IsoStats isoStats;
	if( sampleData )
	{
		SparseNodeData< ProjectiveData< TotalPointSampleData , Real > , IsotropicUIntPack< Dim , DataSig > > _sampleData = tree.template setDataField< DataSig , false >( *samples , *sampleData , (DensityEstimator*)NULL );
		for( const RegularTreeNode< Dim , FEMTreeNodeData >* n = tree.tree().nextNode() ; n ; n=tree.tree().nextNode( n ) )
		{
			ProjectiveData< TotalPointSampleData , Real >* clr = _sampleData( n );
			if( clr ) (*clr) *= (Real)pow( DataX.value , tree.depth( n ) );
		}
		isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< TotalPointSampleData >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , &_sampleData , solution , isoValue , *mesh , SetVertex , !LinearFit.set , !NonManifold.set , PolygonMesh.set , false );
	}
#if defined( __GNUC__ ) && __GNUC__ < 5
	#warning "you've got me gcc version<5"
	else isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< TotalPointSampleData >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , (SparseNodeData< ProjectiveData< TotalPointSampleData , Real > , IsotropicUIntPack< Dim , DataSig > > *)NULL , solution , isoValue , *mesh , SetVertex , !LinearFit.set , !NonManifold.set , PolygonMesh.set , false );
#else // !__GNUC__ || __GNUC__ >=5
	else isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< TotalPointSampleData >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , NULL , solution , isoValue , *mesh , SetVertex , !LinearFit.set , !NonManifold.set , PolygonMesh.set , false );
#endif // __GNUC__ || __GNUC__ < 4
	messageWriter( "Vertices / Polygons: %d / %d\n" , mesh->outOfCorePointCount()+mesh->inCorePoints.size() , mesh->polygonCount() );
	//messageWriter( "Corners / Vertices / Edges / Surface / Set Table / Copy Finer: %.1f / %.1f / %.1f / %.1f / %.1f / %.1f (s)\n" , isoStats.cornersTime , isoStats.verticesTime , isoStats.edgesTime , isoStats.surfaceTime , isoStats.setTableTime , isoStats.copyFinerTime );
	if( PolygonMesh.set ) profiler.dumpOutput2( comments , "#         Got polygons:" );
	else                  profiler.dumpOutput2( comments , "#        Got triangles:" );

	// to do mesh
	if(mesh->polygonCount() > 0)
	{
		mesh->resetIterator();
		my_io.proc_buffers.num_pts = mesh->outOfCorePointCount() + (int)mesh->inCorePoints.size();
		my_io.proc_buffers.num_primitives = (unsigned int)mesh->polygonCount();
		my_io.proc_buffers.pos_pt = new __Coord[my_io.proc_buffers.num_pts];
		if(Normals.set)
			my_io.proc_buffers.nrl_pt = new __Coord[my_io.proc_buffers.num_pts];
		my_io.proc_buffers.index_buffer = new unsigned int[my_io.proc_buffers.num_primitives * Dim];

		typename Vertex::Transform _inverse_xForm(iXForm);
		unsigned int vtxcnt = 0;
		for (int i = 0; i < (int)mesh->inCorePoints.size(); i++, vtxcnt++)
		{
			Vertex vertex = _inverse_xForm(mesh->inCorePoints[i]);

			float* pos = (float*)&my_io.proc_buffers.pos_pt[i];
			float* nrl = Normals.set? (float*)&my_io.proc_buffers.nrl_pt[i] : NULL;
			SetOutBuf(pos, nrl, vertex);
		}
		for (int i = 0; i < mesh->outOfCorePointCount(); i++, vtxcnt++)
		{
			Vertex pt;
			mesh->nextOutOfCorePoint(pt);

			Vertex vertex = _inverse_xForm(pt);

			float* pos = (float*)&my_io.proc_buffers.pos_pt[vtxcnt];
			float* nrl = Normals.set ? (float*)&my_io.proc_buffers.nrl_pt[i] : NULL;
			SetOutBuf(pos, nrl, vertex);
		}

		unsigned int triCnt = 0;
		std::vector< CoredVertexIndex > polygon;
		for (int i = 0; i < (int)my_io.proc_buffers.num_primitives; i++)
		{
			mesh->nextPolygon(polygon);

			for (int j = 0; j < (int)polygon.size(); j++)
			{
				if (polygon[j].inCore)
					my_io.proc_buffers.index_buffer[i * polygon.size() + j]
						= polygon[j].idx;
				else
					my_io.proc_buffers.index_buffer[i * polygon.size() + j]
						= polygon[j].idx + (int)mesh->inCorePoints.size();
			}
		}  // for, write faces
	}

	std::vector< std::string > noComments;
	if (Out.set)
	{
		if (!PlyWritePolygons< Vertex, Real, Dim >(Out.value, mesh, ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE, NoComments.set ? noComments : comments, iXForm))
			ERROR_OUT("Could not write mesh to: %s", Out.value);
	}

	delete mesh;
}

template< class Real , typename ... SampleData , unsigned int ... FEMSigs , typename __Coord >
void Execute( MyInOutStream<__Coord>& my_inout_stream, UIntPack< FEMSigs ... > )
{
	static const int Dim = sizeof ... ( FEMSigs );
	typedef UIntPack< FEMSigs ... > Sigs;
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > Degrees;
	typedef UIntPack< FEMDegreeAndBType< NORMAL_DEGREE , DerivativeBoundary< FEMSignature< FEMSigs >::BType , 1 >::BType >::Signature ... > NormalSigs;
	static const unsigned int DataSig = FEMDegreeAndBType< DATA_DEGREE , BOUNDARY_FREE >::Signature;
	typedef typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
	typedef typename FEMTree< Dim , Real >::template InterpolationInfo< Real , 0 > InterpolationInfo;
	typedef PointStreamNormal< Real , Dim > NormalPointSampleData;
	typedef MultiPointStreamData< Real , SampleData ... > AdditionalPointSampleData;
	typedef MultiPointStreamData< Real , NormalPointSampleData , AdditionalPointSampleData > TotalPointSampleData;
	typedef InputPointStreamWithData< Real , Dim , TotalPointSampleData > InputPointStream;
	typedef TransformedInputPointStreamWithData< Real , Dim , TotalPointSampleData > XInputPointStream;
	std::vector< std::string > comments;
	messageWriter( comments , "*************************************************************\n" );
	messageWriter( comments , "*************************************************************\n" );
	messageWriter( comments , "** Running Screened Poisson Reconstruction (Version %s) **\n" , VERSION );
	messageWriter( comments , "*************************************************************\n" );
	messageWriter( comments , "*************************************************************\n" );

	XForm< Real , Dim+1 > xForm , iXForm;
	if( Transform.set )
	{
		FILE* fp = fopen( Transform.value , "r" );
		if( !fp )
		{
			WARN( "Could not read x-form from: %s" , Transform.value );
			xForm = XForm< Real , Dim+1 >::Identity();
		}
		else
		{
			for( int i=0 ; i<Dim+1 ; i++ ) for( int j=0 ; j<Dim+1 ; j++ )
			{
				float f;
				if( fscanf( fp , " %f " , &f )!=1 ) ERROR_OUT( "Failed to read xform" );
				xForm(i,j) = (Real)f;
			}
			fclose( fp );
		}
	}
	else xForm = XForm< Real , Dim+1 >::Identity();

	char str[1024];
	for( int i=0 ; params[i] ; i++ )
		if( params[i]->set )
		{
			params[i]->writeValue( str );
			if( strlen( str ) ) messageWriter( comments , "\t--%s %s\n" , params[i]->name , str );
			else                messageWriter( comments , "\t--%s\n" , params[i]->name );
		}

	double startTime = Time();
	Real isoValue = 0;

	FEMTree< Dim , Real > tree( MEMORY_ALLOCATOR_BLOCK_SIZE );
	tree.aux_vfield_mode = Confidence.value > 0 &&
		my_inout_stream.nrl_aux_pts.size() > 0;

	FEMTreeProfiler< Dim , Real > profiler( tree );

	if( Depth.set && Width.value>0 )
	{
		WARN( "Both --%s and --%s set, ignoring --%s" , Depth.name , Width.name , Width.name );
		Width.value = 0;
	}

	int pointCount;

	Real pointWeightSum;
	std::vector< typename FEMTree< Dim , Real >::PointSample >* samples = new std::vector< typename FEMTree< Dim , Real >::PointSample >();
	std::vector< TotalPointSampleData >* sampleData = NULL;
	DensityEstimator* density = NULL;
	SparseNodeData< Point< Real , Dim > , NormalSigs >* normalInfo = NULL;
	Real targetValue = (Real)0.5; // implicit function 에서 iso-value 의 기준 (이것을 기준으로 level set surface 정의)

	// Read in the samples (and (optional) color data or confidence data)
	{
		profiler.start();
		InputPointStream* pointStream;
		sampleData = new std::vector< TotalPointSampleData >();
		// inCore case, pointStream refs inCorePoints
		std::vector< std::pair< Point< Real, Dim >, TotalPointSampleData > > inCorePoints;

		if (my_inout_stream.pos_pts.size() + my_inout_stream.pos_aux_pts.size() == 0 && In.set)
		{
			char* ext = GetFileExtension(In.value);
			if (InCore.set)
			{
				InputPointStream *_pointStream;
				if (!strcasecmp(ext, "bnpts")) 
					_pointStream = new BinaryInputPointStreamWithData< Real, Dim, TotalPointSampleData >(In.value, TotalPointSampleData::ReadBinary);
				else if (!strcasecmp(ext, "ply")) 
					_pointStream = new PLYInputPointStreamWithData< Real, Dim, TotalPointSampleData >(In.value, TotalPointSampleData::PlyReadProperties(), TotalPointSampleData::PlyReadNum, TotalPointSampleData::ValidPlyReadProperties);
				else                                    
					_pointStream = new ASCIIInputPointStreamWithData< Real, Dim, TotalPointSampleData >(In.value, TotalPointSampleData::ReadASCII);
				
				Point< Real, Dim > p;
				TotalPointSampleData d;
				while (_pointStream->nextPoint(p, d)) 
					inCorePoints.push_back(std::pair< Point< Real, Dim >, TotalPointSampleData >(p, d));
				
				delete _pointStream;

				pointStream = new MemoryInputPointStreamWithData< Real, Dim, TotalPointSampleData >(inCorePoints.size(), &inCorePoints[0]);
			}
			else
			{
				if (!strcasecmp(ext, "bnpts")) pointStream = new BinaryInputPointStreamWithData< Real, Dim, TotalPointSampleData >(In.value, TotalPointSampleData::ReadBinary);
				else if (!strcasecmp(ext, "ply")) pointStream = new    PLYInputPointStreamWithData< Real, Dim, TotalPointSampleData >(In.value, TotalPointSampleData::PlyReadProperties(), TotalPointSampleData::PlyReadNum, TotalPointSampleData::ValidPlyReadProperties);
				else                                    pointStream = new  ASCIIInputPointStreamWithData< Real, Dim, TotalPointSampleData >(In.value, TotalPointSampleData::ReadASCII);
			}

			delete[] ext;
		}
		else
		{
			for (size_t i = 0; i < my_inout_stream.pos_pts.size(); i++)
			{
				float* pos_coord = (float*)&my_inout_stream.pos_pts[i];
				float* nrl_coord = (float*)&my_inout_stream.nrl_pts[i];
				Point< Real, Dim > p, n;
				for (int ele = 0; ele < Dim; ele++)
				{
					p.coords[ele] = Real(pos_coord[ele]);
					n.coords[ele] = Real(nrl_coord[ele]);
				}

				NormalPointSampleData n_data;
				n_data.data = n;

				TotalPointSampleData d;
				std::get<0>(d.data) = n_data;

				inCorePoints.push_back(std::pair< Point< Real, Dim >, TotalPointSampleData >(p, d));
			}

			for (size_t i = 0; i < my_inout_stream.pos_aux_pts.size(); i++)
			{
				float* pos_coord = (float*)&my_inout_stream.pos_aux_pts[i];
				float* nrl_coord = (float*)&my_inout_stream.nrl_aux_pts[i];
				Point< Real, Dim > p, n;
				for (int ele = 0; ele < Dim; ele++)
				{
					p.coords[ele] = Real(pos_coord[ele]);
					n.coords[ele] = Real(nrl_coord[ele]);
				}

				NormalPointSampleData n_data;
				n_data.data = n;

				TotalPointSampleData d;
				std::get<0>(d.data) = n_data;

				inCorePoints.push_back(std::pair< Point< Real, Dim >, TotalPointSampleData >(p, d));
			}

			pointStream = new MemoryInputPointStreamWithData< Real, Dim, TotalPointSampleData >(inCorePoints.size(), &inCorePoints[0]);
		}

		typename TotalPointSampleData::Transform _xForm(xForm);
		XInputPointStream _pointStream([&](Point< Real, Dim >& p, TotalPointSampleData& d) { p = xForm * p, d = _xForm(d); }, *pointStream);
		if (Width.value > 0)
			xForm = GetPointXForm< Real, Dim >(_pointStream, Width.value, (Real)(Scale.value > 0 ? Scale.value : 1.), Depth.value) * xForm;
		else
			xForm = Scale.value > 0 ? GetPointXForm< Real, Dim >(_pointStream, (Real)Scale.value) * xForm : xForm;

		// make a tree structure of sample nodes
		{
			typename TotalPointSampleData::Transform _xForm(xForm);
			XInputPointStream _pointStream([&](Point< Real, Dim >& p, TotalPointSampleData& d) { p = xForm * p, d = _xForm(d); }, *pointStream);

			auto ProcessDataWithConfidence = [&](const Point< Real, Dim >& p, TotalPointSampleData& d)
			{
				// std::get< 0 >(d.data).data means normal ... Point< Real, Dim > n
				Real l = (Real)Length(std::get< 0 >(d.data).data);
				if (!l || l != l) return (Real)-1.; // l != l means NAN
				return (Real)pow(l, Confidence.value);
			};
			auto ProcessData = [](const Point< Real, Dim >& p, TotalPointSampleData& d)
			{
				Real l = (Real)Length(std::get< 0 >(d.data).data);
				if (!l || l != l) return (Real)-1.;
				std::get< 0 >(d.data).data /= l; // normalized
				return (Real)1.;
			};
			// modified by dojo 
			// Samples and Tree Node 생성 중 호출
			auto ProcessDataWithAuxConfidence = [&](const Point< Real, Dim >& p, const int idx, bool& is_confidence_point, 
				TotalPointSampleData& d)
			{
				// std::get< 0 >(d.data).data means normal ... Point< Real, Dim > n
				Real len = (Real)Length(std::get< 0 >(d.data).data);
				if (!len || len != len) return (Real)-1.; // l != l means NAN

				if (idx >= my_inout_stream.pos_pts.size())
				{
					is_confidence_point = false;
					return (Real)pow(len, Confidence.value);
				}

				is_confidence_point = true;
				std::get< 0 >(d.data).data /= len; // normalized
				return (Real)1.;
			};
			if (Confidence.value > 0)
			{
				my_inout_stream.pos_aux_pts.size() > 0 ?
					pointCount = FEMTreeInitializer< Dim, Real >::template Initialize< TotalPointSampleData >(tree.spaceRoot(), _pointStream,
						Depth.value, *samples, *sampleData, true, tree.nodeAllocator, tree.initializer(),
						ProcessDataWithAuxConfidence) :
					pointCount = FEMTreeInitializer< Dim, Real >::template Initialize< TotalPointSampleData >(tree.spaceRoot(), _pointStream,
						Depth.value, *samples, *sampleData, true, tree.nodeAllocator, tree.initializer(),
						ProcessDataWithConfidence);
			}
			else
				pointCount = FEMTreeInitializer< Dim, Real >::template Initialize< TotalPointSampleData >(tree.spaceRoot(), _pointStream, 
					Depth.value, *samples, *sampleData, true, tree.nodeAllocator, tree.initializer(), 
					ProcessData);
		}
		iXForm = xForm.inverse();

		delete pointStream;

		messageWriter("Input Points / Samples: %d / %d\n", pointCount, samples->size());
		profiler.dumpOutput2(comments, "# Read input into tree:");
	}

	int kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value; // to do //
	if( kernelDepth>Depth.value )
	{
		WARN( "%s can't be greater than %s: %d <= %d" , KernelDepth.name , Depth.name , KernelDepth.value , Depth.value );
		kernelDepth = Depth.value;
	}

	DenseNodeData< Real , Sigs > solution;
	{
		DenseNodeData< Real , Sigs > constraints;
		InterpolationInfo* iInfo = NULL;
		int solveDepth = Depth.value;

		tree.resetNodeIndices();

		// Get the kernel density estimator
		// 여기에서 (겸시겸시) node 에 대한 sample data 가 처리됨.
		{
			profiler.start();
			// WEIGHT_DEGREE ==> Degree ???
			// weight 를 node 에 splatting 
			density = tree.template setDensityEstimator< WEIGHT_DEGREE >( *samples , kernelDepth , SamplesPerNode.value , 1 );
			profiler.dumpOutput2( comments , "#   Got kernel density:" );
		}

		// Transform the Hermite samples into a vector field
		{
			profiler.start(); 
			normalInfo = new SparseNodeData< Point< Real , Dim > , NormalSigs >();
			// max depth node 에서만 
			// node 에서 splatting 된 weight 를 기반으로 B-spline 으로 normal 이  splating 됨.
			if( ConfidenceBias.value>0 ) 
				*normalInfo = tree.setNormalField( NormalSigs() , *samples , *sampleData , density , pointWeightSum , [&]( Real conf ){ return (Real)( log( conf ) * ConfidenceBias.value / log( 1<<(Dim-1) ) ); } );
			else                         
				*normalInfo = tree.setNormalField( NormalSigs() , *samples , *sampleData , density , pointWeightSum );
			
#pragma omp parallel for
			for( int i=0 ; i<normalInfo->size() ; i++ ) (*normalInfo)[i] *= (Real)-1.;
			profiler.dumpOutput2( comments , "#     Got normal field:" );
			messageWriter( "Point weight / Estimated Area: %g / %g\n" , pointWeightSum , pointCount*pointWeightSum );
		}

		if( !Density.set ) delete density , density = NULL;
		if( DataX.value<=0 || ( !Colors.set && !Normals.set ) ) delete sampleData , sampleData = NULL;

		// Trim the tree and prepare for multigrid
		{
			profiler.start();
			constexpr int MAX_DEGREE = NORMAL_DEGREE > Degrees::Max() ? NORMAL_DEGREE : Degrees::Max();
			tree.template finalizeForMultigrid< MAX_DEGREE >( FullDepth.value , typename FEMTree< Dim , Real >::template HasNormalDataFunctor< NormalSigs >( *normalInfo ) , normalInfo , density );
			profiler.dumpOutput2( comments , "#       Finalized tree:" );
		}
		// Add the FEM constraints
		{
			profiler.start();
			constraints = tree.initDenseNodeData( Sigs() );
			typename FEMIntegrator::template Constraint< Sigs , IsotropicUIntPack< Dim , 1 > , NormalSigs , IsotropicUIntPack< Dim , 0 > , Dim > F;
			unsigned int derivatives2[Dim];
			for( int d=0 ; d<Dim ; d++ ) derivatives2[d] = 0;
			typedef IsotropicUIntPack< Dim , 1 > Derivatives1;
			typedef IsotropicUIntPack< Dim , 0 > Derivatives2;
			for( int d=0 ; d<Dim ; d++ )
			{
				unsigned int derivatives1[Dim];
				for( int dd=0 ; dd<Dim ; dd++ ) derivatives1[dd] = dd==d ?  1 : 0;
				F.weights[d][ TensorDerivatives< Derivatives1 >::Index( derivatives1 ) ][ TensorDerivatives< Derivatives2 >::Index( derivatives2 ) ] = 1;
			}
			tree.addFEMConstraints( F , *normalInfo , constraints , solveDepth );
			profiler.dumpOutput2( comments , "#  Set FEM constraints:" );
		}

		// Free up the normal info
		delete normalInfo , normalInfo = NULL;

		// Add the interpolation constraints
		if( PointWeight.value>0 )
		{
			profiler.start();
			if( ExactInterpolation.set ) 
				iInfo = FEMTree< Dim , Real >::template InitializeExactPointInterpolationInfo< Real , 0 > ( tree , *samples , ConstraintDual< Dim , Real >( targetValue , (Real)PointWeight.value * pointWeightSum ) , SystemDual< Dim , Real >( (Real)PointWeight.value * pointWeightSum ) , true , false );
			else                         
				iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointInterpolationInfo< Real , 0 > ( tree , *samples , ConstraintDual< Dim , Real >( targetValue , (Real)PointWeight.value * pointWeightSum ) , SystemDual< Dim , Real >( (Real)PointWeight.value * pointWeightSum ) , true , 1 );
			
			tree.addInterpolationConstraints( constraints , solveDepth , *iInfo );
			profiler.dumpOutput2( comments , "#Set point constraints:" );
		}

		messageWriter( "Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n" , (int)tree.leaves() , (int)tree.nodes() , (int)tree.ghostNodes() );
		messageWriter( "Memory Usage: %.3f MB\n" , float( MemoryInfo::Usage())/(1<<20) );
		
		// Solve the linear system
		{
			profiler.start();
			typename FEMTree< Dim , Real >::SolverInfo sInfo;
			sInfo.cgDepth = CGDepth.value, sInfo.cascadic = true , sInfo.vCycles = 1 , sInfo.iters = Iters.value , sInfo.cgAccuracy = CGSolverAccuracy.value , sInfo.verbose = Verbose.set , sInfo.showResidual = ShowResidual.set , sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE , sInfo.sliceBlockSize = 1;
			sInfo.baseDepth = BaseDepth.value , sInfo.baseVCycles = BaseVCycles.value;
			typename FEMIntegrator::template System< Sigs , IsotropicUIntPack< Dim , 1 > > F( { 0. , 1. } );
			solution = tree.solveSystem( Sigs() , F , constraints , solveDepth , sInfo , iInfo );
			profiler.dumpOutput2( comments , "# Linear system solved:" );
			if( iInfo ) delete iInfo , iInfo = NULL;
		}
	}

	{
		profiler.start();
		double valueSum = 0 , weightSum = 0;
		typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< Sigs , 0 > evaluator( &tree , solution );
#pragma omp parallel for reduction( + : valueSum , weightSum )
		for( int j=0 ; j<samples->size() ; j++ )
		{
			ProjectiveData< Point< Real , Dim > , Real >& sample = (*samples)[j].sample;
			Real w = sample.weight;
			if (w > 0 && (*samples)[j].node->confidence_flag != CONFI_VF_AUX_ONLY)
			{
				weightSum += w;
				valueSum += evaluator.values(sample.data / sample.weight, omp_get_thread_num(), (*samples)[j].node)[0] * w;
			}
		}
		isoValue = weightSum > 0 ? (Real)( valueSum / weightSum ) : 0;
		//isoValue = 0;
		if( DataX.value<=0 || ( !Colors.set && !Normals.set ) ) delete samples , samples = NULL;
		profiler.dumpOutput( "Got average:" );
		messageWriter( "Iso-Value: %e = %g / %g\n" , isoValue , valueSum , weightSum );
	}
	if( Tree.set )
	{
		FILE* fp = fopen( Tree.value , "wb" );
		if( !fp ) ERROR_OUT( "Failed to open file for writing: %s" , Tree.value );
		FEMTree< Dim , Real >::WriteParameter( fp );
		DenseNodeData< Real , Sigs >::WriteSignatures( fp );
		tree.write( fp );
		solution.write( fp );
		fclose( fp );
	}

	printf("=================> 1\n");
	int res = 0;
	if (LevelSetMap.set)
	{
		Pointer2(Real) values2d = tree.template regularGridEvaluate_2darray< true >(solution, res, -1, PrimalVoxel.set);
		my_inout_stream.proc_buffers.level_set = values2d;
		my_inout_stream.proc_buffers.res = res;
		my_inout_stream.proc_buffers.iso_value = (float)isoValue;
	}
	printf("=================> 2\n");

	Pointer(Real) values = NULL;
	if( VoxelGrid.set )
	{
		if(values == NULL)
			values = tree.template regularGridEvaluate< true >(solution, res, -1, PrimalVoxel.set);

		FILE* fp = fopen( VoxelGrid.value , "wb" );
		if (!fp)
		{
			WARN("Failed to open voxel file for writing: %s", VoxelGrid.value);
		}
		else
		{
			profiler.start();
			int num_voxels = 1;
			for (int i = 0; i < Dim; i++) num_voxels *= res;

#pragma omp parallel for
			for( int i=0 ; i< num_voxels; i++ ) values[i] -= isoValue;
			profiler.dumpOutput( "Got voxel grid:" );
			fwrite( &res , sizeof(int) , 1 , fp );
			if( typeid(Real)==typeid(float) ) fwrite( values , sizeof(float) , num_voxels, fp );
			else
			{
				float *fValues = new float[num_voxels];
				for( int i=0 ; i< num_voxels; i++ ) fValues[i] = float( values[i] );
				fwrite( fValues , sizeof(float) , num_voxels, fp );
				delete[] fValues;
			}
			fclose( fp );
		}
	}

	if (!SkipContouring.set)
	{
#define __EXTRACT_MESH
#ifdef __EXTRACT_MESH
		if (Dim == 3) // Out.set
		{
			if (Normals.set)
			{
				if (Density.set)
				{
					typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, PointStreamValue< Real >, AdditionalPointSampleData > > Vertex;
					std::function< void(Vertex&, Point< Real, Dim >, Real, TotalPointSampleData) > SetVertex = [](Vertex& v, Point< Real, Dim > p, Real w, TotalPointSampleData d) { v.point = p, std::get< 0 >(v.data.data) = std::get< 0 >(d.data), std::get< 1 >(v.data.data).data = w, std::get< 2 >(v.data.data) = std::get< 1 >(d.data); };
					std::function< void(float*, float*, Vertex&) > SetOutBuf = [](float* pos, float* nrl, Vertex& v)
					{
						for (int ele = 0; ele < Dim; ele++)
						{
							pos[ele] = (float)v.point[ele];
							nrl[ele] = std::get<0>(v.data.data).data[ele];
						}
					};
					ExtractMesh< Dim, Vertex >(my_inout_stream, UIntPack< FEMSigs ... >(), std::tuple< SampleData ... >(), tree, solution, isoValue, samples, sampleData, density, SetVertex, comments, iXForm, SetOutBuf);
				}
				else
				{
					typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, AdditionalPointSampleData > > Vertex;
					std::function< void(Vertex&, Point< Real, Dim >, Real, TotalPointSampleData) > SetVertex = [](Vertex& v, Point< Real, Dim > p, Real w, TotalPointSampleData d) { v.point = p, std::get< 0 >(v.data.data) = std::get< 0 >(d.data), std::get< 1 >(v.data.data) = std::get< 1 >(d.data);  };
					std::function< void(float*, float*, Vertex&) > SetOutBuf = [](float* pos, float* nrl, Vertex& v)
					{
						for (int ele = 0; ele < Dim; ele++)
						{
							pos[ele] = (float)v.point[ele];
							nrl[ele] = std::get<0>(v.data.data).data[ele];
						}
					};
					ExtractMesh< Dim, Vertex >(my_inout_stream, UIntPack< FEMSigs ... >(), std::tuple< SampleData ... >(), tree, solution, isoValue, samples, sampleData, density, SetVertex, comments, iXForm, SetOutBuf);
				}
			}
			else
			{
				if (Density.set)
				{
					typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamValue< Real >, AdditionalPointSampleData > > Vertex;
					std::function< void(Vertex&, Point< Real, Dim >, Real, TotalPointSampleData) > SetVertex = [](Vertex& v, Point< Real, Dim > p, Real w, TotalPointSampleData d) { v.point = p, std::get< 0 >(v.data.data).data = w, std::get< 1 >(v.data.data) = std::get< 1 >(d.data); };
					std::function< void(float*, float*, Vertex&) > SetOutBuf = [](float* pos, float* nrl, Vertex& v)
					{
						for (int ele = 0; ele < Dim; ele++)
						{
							pos[ele] = (float)v.point[ele];
						}
					};
					ExtractMesh< Dim, Vertex >(my_inout_stream, UIntPack< FEMSigs ... >(), std::tuple< SampleData ... >(), tree, solution, isoValue, samples, sampleData, density, SetVertex, comments, iXForm, SetOutBuf);
				}
				else
				{
					typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, AdditionalPointSampleData > > Vertex;
					std::function< void(Vertex&, Point< Real, Dim >, Real, TotalPointSampleData) > SetVertex = [](Vertex& v, Point< Real, Dim > p, Real w, TotalPointSampleData d) { v.point = p, std::get< 0 >(v.data.data) = std::get< 1 >(d.data); };
					std::function< void(float*, float*, Vertex&) > SetOutBuf = [](float* pos, float* nrl, Vertex& v)
					{
						for (int ele = 0; ele < Dim; ele++)
						{
							pos[ele] = (float)v.point[ele];
						}
					};
					ExtractMesh< Dim, Vertex >(my_inout_stream, UIntPack< FEMSigs ... >(), std::tuple< SampleData ... >(), tree, solution, isoValue, samples, sampleData, density, SetVertex, comments, iXForm, SetOutBuf);
				}
			}
		}
		else if (Dim == 2)
#endif
		{
			if (values == NULL)
				values = tree.template regularGridEvaluate< true >(solution, res, -1, PrimalVoxel.set);

			std::vector<__float2> lines_list;
			extract_iso_contour<Real>(lines_list, (float)isoValue, values, res, res);

			my_inout_stream.proc_buffers.index_buffer = new unsigned int[lines_list.size()];
			my_inout_stream.proc_buffers.pos_pt = new __Coord[lines_list.size()];
			my_inout_stream.proc_buffers.num_pts = (unsigned int)lines_list.size();
			my_inout_stream.proc_buffers.num_primitives = (unsigned int)lines_list.size() / 2;

			for (size_t i = 0; i < lines_list.size(); i++)
			{
				Point< Real, Dim > p;
				p.coords[0] = (Real)lines_list[i].x / (Real)res;
				p.coords[1] = (Real)lines_list[i].y / (Real)res;
				p = iXForm * p;

				float* pos = (float*)&my_inout_stream.proc_buffers.pos_pt[i];
				for (int j = 0; j < Dim; j++)
					pos[j] = (float)p.coords[j];
				my_inout_stream.proc_buffers.index_buffer[i] = (unsigned int)i;
			}
		}
	}

	printf("=================> 3\n");
	//if (!LevelSetMap.set) 
		DeletePointer(values);
	printf("=================> 4\n");

	if (sampleData) { delete sampleData; sampleData = NULL; }
	if( density ) delete density , density = NULL;
	messageWriter( comments , "#          Total Solve: %9.1f (s), %9.1f (MB)\n" , Time()-startTime , FEMTree< Dim , Real >::MaxMemoryUsage() );
}

#ifndef FAST_COMPILE
template< unsigned int Dim , class Real , typename ... SampleData , typename __Coord >
void Execute(MyInOutStream<__Coord>& my_inout_stream)
{
	switch( BType.value )
	{
		case BOUNDARY_FREE+1:
		{
			switch( Degree.value )
			{
				case 1: return Execute< Real , SampleData ... >( my_inout_stream , IsotropicUIntPack< Dim , FEMDegreeAndBType< 1 , BOUNDARY_FREE >::Signature >() );
				case 2: return Execute< Real , SampleData ... >( my_inout_stream , IsotropicUIntPack< Dim , FEMDegreeAndBType< 2 , BOUNDARY_FREE >::Signature >() );
//				case 3: return Execute< Real , SampleData ... >( my_input_stream , IsotropicUIntPack< Dim , FEMDegreeAndBType< 3 , BOUNDARY_FREE >::Signature >() );
//				case 4: return Execute< Real , SampleData ... >( my_input_stream , IsotropicUIntPack< Dim , FEMDegreeAndBType< 4 , BOUNDARY_FREE >::Signature >() );
				default: ERROR_OUT( "Only B-Splines of degree 1 - 2 are supported" );
			}
		}
		case BOUNDARY_NEUMANN+1:
		{
			switch( Degree.value )
			{
				case 1: return Execute< Real , SampleData ... >( my_inout_stream , IsotropicUIntPack< Dim , FEMDegreeAndBType< 1 , BOUNDARY_NEUMANN >::Signature >() );
				case 2: return Execute< Real , SampleData ... >( my_inout_stream , IsotropicUIntPack< Dim , FEMDegreeAndBType< 2 , BOUNDARY_NEUMANN >::Signature >() );
//				case 3: return Execute< Real , SampleData ... >( my_input_stream , IsotropicUIntPack< Dim , FEMDegreeAndBType< 3 , BOUNDARY_NEUMANN >::Signature >() );
//				case 4: return Execute< Real , SampleData ... >( my_input_stream , IsotropicUIntPack< Dim , FEMDegreeAndBType< 4 , BOUNDARY_NEUMANN >::Signature >() );
				default: ERROR_OUT( "Only B-Splines of degree 1 - 2 are supported" );
			}
		}
		case BOUNDARY_DIRICHLET+1:
		{
			switch( Degree.value )
			{
			case 1: return Execute< Real, SampleData ... >(my_inout_stream, IsotropicUIntPack< Dim, FEMDegreeAndBType< 1, BOUNDARY_DIRICHLET >::Signature >());
			case 2: return Execute< Real, SampleData ... >(my_inout_stream, IsotropicUIntPack< Dim, FEMDegreeAndBType< 2, BOUNDARY_DIRICHLET >::Signature >());
//			case 3: return Execute< Real , SampleData ... >( my_input_stream , IsotropicUIntPack< Dim , FEMDegreeAndBType< 3 , BOUNDARY_DIRICHLET >::Signature >() );
//			case 4: return Execute< Real , SampleData ... >( my_input_stream , IsotropicUIntPack< Dim , FEMDegreeAndBType< 4 , BOUNDARY_DIRICHLET >::Signature >() );
			default: ERROR_OUT( "Only B-Splines of degree 1 - 2 are supported" );
			}
		}
		default: ERROR_OUT( "Not a valid boundary type: %d" , BType.value );
	}
}
#endif // !FAST_COMPILE

int main( int argc , char* argv[] )
{
	Timer timer;
#ifdef ARRAY_DEBUG
	WARN( "Array debugging enabled" );
#endif // ARRAY_DEBUG

	//cmdLineParse( argc-1 , &argv[1] , params );
	if( MaxMemoryGB.value>0 ) SetPeakMemoryMB( MaxMemoryGB.value<<10 );

#ifdef _OPENMP
	Threads.value = omp_get_num_threads();
#endif
#ifdef	_DEBUG
	Threads.value = 1;
#endif

	std::string str_in = std::string("D:\\Data\\Experiments 2019\\piston_half_crop\\test_cube.ply");
	std::string str_out = std::string("test_out.ply");
	In.value = (char*)str_in.c_str();
	In.set = true;
	Out.value = (char*)str_out.c_str();
	Out.set = true;
	PointWeight.value = 10;
	PointWeight.set = true;

	Scale.value = 1.1f;
	Scale.set = true;
	
	Depth.value = 5;
	Depth.set = true;
	FullDepth.value = 3; // 이전 버전에서는 
	FullDepth.set = true;
	KernelDepth.value = Depth.value;
	KernelDepth.set = true;
	CGDepth.value = 0;
	CGDepth.set = true;
	BaseDepth.value = FullDepth.value;
	BaseDepth.set = true;

	Verbose.set = true;
	LinearFit.set = true;


	InCore.set = true;



	omp_set_num_threads( Threads.value > 1 ? Threads.value : 1 );
	messageWriter.echoSTDOUT = Verbose.set;

	if( !In.set )
	{
		ShowUsage( argv[0] );
		return 0;
	}
	if( DataX.value<=0 ) Normals.set = Colors.set = false;
	if( BaseDepth.value>FullDepth.value )
	{
		if( BaseDepth.set ) WARN( "Base depth must be smaller than full depth: %d <= %d" , BaseDepth.value , FullDepth.value );
		BaseDepth.value = FullDepth.value;
	}

#ifdef USE_DOUBLE
	typedef double Real;
#else // !USE_DOUBLE
	typedef float  Real;
#endif // USE_DOUBLE

	std::vector<__float3> __dummy;
	MyInOutStream< __float3> my_is(__dummy, __dummy, __dummy, __dummy);

#ifdef FAST_COMPILE
	static const int Degree = DEFAULT_FEM_DEGREE;
	static const BoundaryType BType = DEFAULT_FEM_BOUNDARY;
	typedef IsotropicUIntPack< DIMENSION , FEMDegreeAndBType< Degree , BType >::Signature > FEMSigs;
	WARN( "Compiled for degree-%d, boundary-%s, %s-precision _only_" , Degree , BoundaryNames[ BType ] , sizeof(DefaultFloatType)==4 ? "single" : "double" );
	if( !PointWeight.set ) PointWeight.value = DefaultPointWeightMultiplier*Degree;
	if( Colors.set ) Execute< Real , PointStreamColor< DefaultFloatType > >( my_is , FEMSigs() );
	else             Execute< Real >( my_is , FEMSigs() );
#else // !FAST_COMPILE
	if( !PointWeight.set ) PointWeight.value = DefaultPointWeightMultiplier*Degree.value;
	if( Colors.set ) Execute< DIMENSION , Real , PointStreamColor< float > >( my_is );
	else             Execute< DIMENSION , Real >( my_is );
#endif // FAST_COMPILE
	if( Performance.set )
	{
		printf( "Time (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
		printf( "Peak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );
	}

	In.value = NULL;
	In.set = false;
	Out.value = NULL;
	Out.set = false;
	return EXIT_SUCCESS;
}

bool ScreenedPoissonSurface2D(__ProcBuffers<__float2>* pOut,
	const std::vector<__float2>& pos_pts, const std::vector<__float2>& nrl_pts,
	const std::vector<__float2>& pos_aux_pts, const std::vector<__float2>& nrl_aux_pts,
	std::vector<char*>& cmd_sps_params, std::vector<void*>& additional_params)
{
	cmdLineParse((int)cmd_sps_params.size(), &cmd_sps_params[0], params);
	
	InCore.set = true;
	Verbose.set = true;
	Scale.value = max(Scale.value, 1.1f);

	if (MaxMemoryGB.value > 0) SetPeakMemoryMB(MaxMemoryGB.value << 10);
#ifdef _OPENMP
	Threads.value = omp_get_num_threads();
#endif
#ifdef	_DEBUG
	Threads.value = 1;
#endif

	typedef float  Real;
	// to do
	MyInOutStream<__float2> my_inout_stream(pos_pts, nrl_pts, pos_aux_pts, nrl_aux_pts);
	Execute< 2, Real >(my_inout_stream);

	*pOut = my_inout_stream.proc_buffers;
	my_inout_stream.proc_buffers.pos_pt = NULL;
	my_inout_stream.proc_buffers.nrl_pt = NULL;
	my_inout_stream.proc_buffers.index_buffer = NULL;
	my_inout_stream.proc_buffers.level_set = NULL;

	return true;
}

bool ScreenedPoissonSurface3D(__ProcBuffers< __float3>* pOut,
	const std::vector<__float3>& pos_pts, const std::vector<__float3>& nrl_pts,
	const std::vector<__float3>& pos_aux_pts, const std::vector<__float3>& nrl_aux_pts,
	std::vector<char*>& cmd_sps_params, std::vector<void*>& additional_params)
{
	cmdLineParse((int)cmd_sps_params.size(), &cmd_sps_params[0], params);

	InCore.set = true;
	
	//Scale.value = (float)1.1;
	//Scale.set = true;
	//Confidence.value = 0;
	//Confidence.set = false;
	//PointWeight.value = 10;
	//PointWeight.set = true;
	//Depth.value = 5;
	//Depth.set = true;
	//FullDepth.value = 0; // 이전 버전에서는 
	//FullDepth.set = true;
	//KernelDepth.value = Depth.value;
	//KernelDepth.set = true;
	//CGDepth.value = 6;
	//CGDepth.set = true;
	//BaseDepth.value = 0;
	//BaseDepth.set = true;
	//LinearFit.set = true;
	////Normals.set = true;
	//BType.value = 3;
	//BType.set = true;
	
	Verbose.set = true;
	Scale.value = max(Scale.value, 1.1f);

	omp_set_num_threads(Threads.value > 1 ? Threads.value : 1);
	std::cout << "process cores : " << Threads.value << std::endl;

	if (pos_aux_pts.size() > 0)
	{
		std::cout << "@@@ aux points force to set Density as false" << std::endl;
		Density.set = false;
	}

	if (MaxMemoryGB.value > 0) SetPeakMemoryMB(MaxMemoryGB.value << 10);
#ifdef _OPENMP
	Threads.value = omp_get_num_threads();
#endif
#ifdef	_DEBUG
	Threads.value = 1;
#endif
	omp_set_num_threads(Threads.value > 1 ? Threads.value : 1);
	std::cout << "process cores : " << Threads.value << std::endl;

	MyInOutStream<__float3> my_inout_stream(pos_pts, nrl_pts, pos_aux_pts, nrl_aux_pts);

	typedef float  Real;

	Timer timer;

	if (!PointWeight.set)
		PointWeight.value = DefaultPointWeightMultiplier * Degree.value;

	if (BaseDepth.value > FullDepth.value)
	{
		if (BaseDepth.set) WARN("Base depth must be smaller than full depth: %d <= %d", BaseDepth.value, FullDepth.value);
		BaseDepth.value = FullDepth.value;
	}

	//{
	//	using namespace std;
	//	cout << "conf points : " << pos_pts.size() << endl;
	//	cout << "aux field points : " << pos_aux_pts.size() << endl;
	//
	//	auto gm = [](const __float3& v)
	//	{
	//		return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	//	};
	//	auto mmx = [&gm](float& _min, float& _max, const vector<__float3>& nrls)
	//	{
	//		_min = FLT_MAX, _max = 0;
	//		for (int i = 0; i < (int)nrls.size(); i++)
	//		{
	//			_min = min(_min, gm(nrls[i]));
	//			_max = max(_max, gm(nrls[i]));
	//		}
	//	};
	//	float _min, _max;
	//	mmx(_min, _max, nrl_pts);
	//	cout << "conf point min max : " << _min << ", " << _max << endl;
	//	mmx(_min, _max, nrl_aux_pts);
	//	cout << "aux point min max : " << _min << ", " << _max << endl;
	//}

	Execute< DIMENSION/*3*/, Real >(my_inout_stream);

	//ShowUsage("my test");

	*pOut = my_inout_stream.proc_buffers;
	my_inout_stream.proc_buffers.pos_pt = NULL;
	my_inout_stream.proc_buffers.nrl_pt = NULL;
	my_inout_stream.proc_buffers.index_buffer = NULL;
	my_inout_stream.proc_buffers.level_set = NULL;

	if (Performance.set)
	{
		printf("Time (Wall/CPU): %.2f / %.2f\n", timer.wallTime(), timer.cpuTime());
		printf("Peak Memory (MB): %d\n", MemoryInfo::PeakMemoryUsageMB());
	}

	return true;
}