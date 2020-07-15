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

template< class Real , unsigned int DataDegree , unsigned int ... DataDegrees > typename std::enable_if< sizeof ... ( DataDegrees )==0 >::type __SetBSplineComponentValues( const Real* position , const Real* start , Real width , double* values , unsigned int stride )
{
	Polynomial< DataDegree >::BSplineComponentValues( ( position[0] - start[0] ) / width , values );
}
template< class Real , unsigned int DataDegree , unsigned int ... DataDegrees > typename std::enable_if< sizeof ... ( DataDegrees )!=0 >::type __SetBSplineComponentValues( const Real* position , const Real* start , Real width , double* values , unsigned int stride )
{
	Polynomial< DataDegree >::BSplineComponentValues( ( position[0] - start[0] ) / width , values );
	__SetBSplineComponentValues< Real , DataDegrees ... >( position+1 , start+1 , width , values + stride , stride );
}


// evaluate the result of splatting along a plane and then evaluating at a point on the plane.
template< unsigned int Degree > double GetScaleValue( void )
{
	double centerValues[Degree+1];
	Polynomial< Degree >::BSplineComponentValues( 0.5 , centerValues );
	double scaleValue = 0;
	for( int i=0 ; i<=Degree ; i++ ) scaleValue += centerValues[i] * centerValues[i];
	return 1./ scaleValue;
}
template< unsigned int Dim , class Real >
template< unsigned int WeightDegree >
void FEMTree< Dim , Real >::_addWeightContribution( DensityEstimator< WeightDegree >& densityWeights , FEMTreeNode* node , Point< Real , Dim > position , PointSupportKey< IsotropicUIntPack< Dim , WeightDegree > >& weightKey , Real weight )
{
	// support element 에 대한 각 element 의 contribution 합이 1 이 되도록 scale
	// 이것이 너무 크면 점 중심으로 energy 가 몰림 (점 중심으로 구들이 생성됨)
	// 너무 작으면 점 밖으로 energy 가 나감
	// 적당한 값은... weight 의 weighed average 값을 취함
	// 아마 density estimation 하는 논문에서 이렇게 하라고 설명되어 있을 듯... 나중에 확인..
	static const double ScaleValue = GetScaleValue< WeightDegree >();
	double values[ Dim ][ BSplineSupportSizes< WeightDegree >::SupportSize ];
	// getNeighbors ==> getNeighbors< 'true' > 에서 true 는... sample node 로부터 생성된 node 외의 모든 이웃 node 를 get 
	// 이 과정에서 새로운 node 가 생성될 수 있으며, 해당 node 의 부모는 null 일수도 있음
	// 이런 의미에서 if(IsActiveNode< Dim >( _node )) 은 sample node 로부터 생성된 node 만 보는 것이고
	// 여기에서 사용되는 if(_node) 은 모든 neighbor 를 보는 것임
	typename FEMTreeNode::template Neighbors< IsotropicUIntPack< Dim , BSplineSupportSizes< WeightDegree >::SupportSize > >& neighbors = weightKey.template getNeighbors< true >( node , nodeAllocator , _NodeInitializer( *this ) );

	// modified by dojo
	// splatting 시 주변 sample / node 에 주는 영향에 대한 것으로
	// 여기서 핵심은
	// conf ==> 비 aux field point sample 에 영향 주도록
	// aux field point sample 은 자기 자신에만 영향 주도록
	// 하는 것.

	densityWeights.reserve( nodeCount() );

	Point< Real , Dim > start;
	Real w; // node 의 width 의미
	_startAndWidth( node , start , w );

	// values ==> spline weight values
	// position 이 node 중심점 대비 shift 된 것에 대한 BSpline...
	for( int dim=0 ; dim<Dim ; dim++ ) Polynomial< WeightDegree >::BSplineComponentValues( ( position[dim]-start[dim] ) / w , values[dim] );

	weight *= (Real)ScaleValue;
	double scratch[Dim+1];
	scratch[0] = weight;
	// modified by dojo
	// node 값이 이웃들로 퍼져나가는 방식 (BSpline support 를 기반으로 splatting)
	// aux 는 sample node 에 충분히 매핑되도록 주어짐... 퍼지지 않게...
	// aux 는 density 에 해당 weight 만 주어지도록...
	// conf 는 aux 외에만 퍼지도록...
	WindowLoop< Dim >::Run
	(
		IsotropicUIntPack< Dim , 0 >() , IsotropicUIntPack< Dim , BSplineSupportSizes< WeightDegree >::SupportSize >() ,
		[&]( int d , int i ){ 
			scratch[d+1] = scratch[d] * values[d][i]; 
		} ,
		[&]( FEMTreeNode* _node ){ 
			// 여기서는 모든 neighbor 를 봄 (sample node 로부터 파생된 node 뿐만 아니라 모든 이웃 node)
			// _node->confidence_flag 가 DEFAULT 일 수 있음 (not set)
			// cf. if(IsActiveNode< Dim >( _node ))
			if (_node) 
			{
				//if (node->confidence_flag == CONFI_VF_AUX_ONLY)
				//{
				//	//if (_node == node) // 자기 자신 일 때만 합...?! 아니면 아예 배제?
				//	//	AddAtomic(densityWeights[_node], (Real)scratch[Dim]);
				//	if (_node->confidence_flag != CONFI_VF_VISITED)
				//		AddAtomic(densityWeights[_node], (Real)scratch[Dim]);
				//}
				//// node is CONFI_VF_VISITED
				//// ISUUE TO DO............. 그냥 confidence node 는 다 splatting 되게 하는 것이 맞는 것 같다..??
				//else if (_node->confidence_flag != CONFI_VF_AUX_ONLY)
					AddAtomic(densityWeights[_node], (Real)scratch[Dim]);
			}
		} ,
		neighbors.neighbors()
	);
}

template< unsigned int Dim , class Real >
template< unsigned int WeightDegree , class PointSupportKey >
Real FEMTree< Dim , Real >::_getSamplesPerNode( const DensityEstimator< WeightDegree >& densityWeights , const FEMTreeNode* node , Point< Real , Dim > position , PointSupportKey& weightKey ) const
{
	Real weight = 0;
	typedef typename PointSupportKey::NeighborType Neighbors;
	double values[ Dim ][ BSplineSupportSizes< WeightDegree >::SupportSize ];
	Neighbors neighbors = weightKey.getNeighbors( node );
	Point< Real , Dim > start;
	Real w;
	_startAndWidth( node , start , w );

	for( int dim=0 ; dim<Dim ; dim++ ) Polynomial< WeightDegree >::BSplineComponentValues( ( position[dim]-start[dim] ) / w , values[dim] );
	double scratch[Dim+1];
	scratch[0] = 1;
	WindowLoop< Dim >::Run
	(
		IsotropicUIntPack< Dim , 0 >() , IsotropicUIntPack< Dim , BSplineSupportSizes< WeightDegree >::SupportSize >() ,
		[&]( int d , int i ){ scratch[d+1] = scratch[d] * values[d][i]; } ,
		[&]( typename Neighbors::Window::data_type node ){ if( node ){ const Real* w = densityWeights( node ) ; if( w ) weight += (Real)( scratch[Dim] * (*w) ); } } ,
		neighbors.neighbors()
	);
	return weight;
}
template< unsigned int Dim , class Real >
template< unsigned int WeightDegree , class PointSupportKey >
void FEMTree< Dim , Real >::_getSampleDepthAndWeight( const DensityEstimator< WeightDegree >& densityWeights , const FEMTreeNode* node , Point< Real , Dim > position , PointSupportKey& weightKey , Real& depth , Real& weight ) const
{
	const FEMTreeNode* temp = node;
	while( _localDepth( temp )>densityWeights.kernelDepth() ) temp = temp->parent;
	weight = _getSamplesPerNode( densityWeights , temp , position , weightKey );
	if( weight>=(Real)1. ) depth = Real( _localDepth( temp ) + log( weight ) / log(double(1<<( Dim-densityWeights.coDimension() ))) );
	else
	{
		Real oldWeight , newWeight;
		oldWeight = newWeight = weight;
		while( newWeight<(Real)1. && temp->parent )
		{
			temp=temp->parent;
			oldWeight = newWeight;
			newWeight = _getSamplesPerNode( densityWeights , temp , position , weightKey );
		}
		depth = Real( _localDepth( temp ) + log( newWeight ) / log( newWeight / oldWeight ) );
	}
	weight = Real( pow( double(1<<( Dim-densityWeights.coDimension() )) , -double(depth) ) );
}
template< unsigned int Dim , class Real >
template< unsigned int WeightDegree , class PointSupportKey >
void FEMTree< Dim , Real >::_getSampleDepthAndWeight( const DensityEstimator< WeightDegree >& densityWeights , Point< Real , Dim > position , PointSupportKey& weightKey , Real& depth , Real& weight ) const
{
	FEMTreeNode* temp;
	Point< Real,  Dim > myCenter;
	for( int d=0 ; d<Dim ; d++ ) myCenter[d] = (Real)0.5;
	Real myWidth = Real( 1. );

	// Get the finest node with depth less than or equal to the splat depth that contains the point
	temp = _spaceRoot;
	while( _localDepth( temp )<densityWeights.kernelDepth() )
	{
		if( !IsActiveNode< Dim >( temp->children ) ) break; // ERROR_OUT( "" );
		int cIndex = FEMTreeNode::ChildIndex( myCenter , position );
		temp = temp->children + cIndex;
		myWidth /= 2;
		for( int d=0 ; d<Dim ; d++ )
			if( (cIndex>>d) & 1 ) myCenter[d] += myWidth/2;
			else                  myCenter[d] -= myWidth/2;
	}
	return _getSampleDepthAndWeight( densityWeights , temp , position , weightKey , depth , weight );
}

template< unsigned int Dim , class Real >
template< bool CreateNodes , class V , unsigned int ... DataSigs >
void FEMTree< Dim , Real >::_splatPointData( FEMTreeNode* node , Point< Real , Dim > position , V v , SparseNodeData< V , UIntPack< DataSigs ... > >& dataInfo , PointSupportKey< UIntPack< FEMSignature< DataSigs >::Degree ... > >& dataKey )
{
	// density 가 있을 경우 density 를 parameter 로 갖는 _splatPointData 에서 호출
	typedef UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportSize ... > SupportSizes;
	double values[ Dim ][ SupportSizes::Max() ];
	typename FEMTreeNode::template Neighbors< UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportSize ... > >& neighbors = dataKey.template getNeighbors< CreateNodes >( node , nodeAllocator , _NodeInitializer( *this ) );

	Point< Real , Dim > start;
	Real w;
	_startAndWidth( node , start , w );

	__SetBSplineComponentValues< Real , FEMSignature< DataSigs >::Degree ... >( &position[0] , &start[0] , w , &values[0][0] , SupportSizes::Max() );
	double scratch[Dim+1];
	scratch[0] = 1;
	// modified by dojo
	// _addWeightContribution 와 동일한 정책
	// aux normal 은 해당 node 에만 set
	WindowLoop< Dim >::Run
	(
		ZeroUIntPack< Dim >() , UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportSize ... >() ,
		[&]( int d , int i ){ 
			scratch[d+1] = scratch[d] * values[d][i]; 
		} ,
		[&]( FEMTreeNode* _node ){ 
			// IsActiveNode ==> 자신을 포함 모든 상층 (부모.. 부모의 부모...) node 가 모두 null 이 아닐 때 true
			// 즉, sample node 로부터 생성된 node 만 취급
			if (IsActiveNode< Dim >(_node))
			{
				//if (node->confidence_flag == CONFI_VF_AUX_ONLY)
				//{
				//	//if (node == _node)
				//	//	AddAtomic(dataInfo[_node], v * (Real)scratch[Dim]);
				//	if (_node->confidence_flag != CONFI_VF_VISITED)
				//		AddAtomic(dataInfo[_node], v * (Real)scratch[Dim]);
				//}
				//else if(_node->confidence_flag != CONFI_VF_AUX_ONLY)
				{
					AddAtomic(dataInfo[_node], v * (Real)scratch[Dim]);
				}
			}
		} ,
		neighbors.neighbors()
	);
}
template< unsigned int Dim , class Real >
template< bool CreateNodes , unsigned int WeightDegree , class V , unsigned int ... DataSigs >
Real FEMTree< Dim , Real >::_splatPointData( const DensityEstimator< WeightDegree >& densityWeights , Point< Real , Dim > position , V v , SparseNodeData< V , UIntPack< DataSigs ... > >& dataInfo , PointSupportKey< IsotropicUIntPack< Dim , WeightDegree > >& weightKey , PointSupportKey< UIntPack< FEMSignature< DataSigs >::Degree ... > >& dataKey , LocalDepth minDepth , LocalDepth maxDepth , int dim , Real depthBias )
{
	// density 가 있을 경우 이것이 호출되며, 현 버전에서는 일단 이것이 호출되어 수행됨
	double dx;
	V _v;
	FEMTreeNode* temp;
	int cnt=0;
	double width;
	Point< Real , Dim > myCenter;
	for( int d=0 ; d<Dim ; d++ ) myCenter[d] = (Real)0.5;
	Real myWidth = (Real)1.;

	temp = _spaceRoot;
	while( _localDepth( temp )<densityWeights.kernelDepth() )
	{
		if( !IsActiveNode< Dim >( temp->children ) ) break;
		int cIndex = FEMTreeNode::ChildIndex( myCenter , position );
		temp = temp->children + cIndex;
		myWidth /= 2;
		for( int d=0 ; d<Dim ; d++ )
			if( (cIndex>>d) & 1 ) myCenter[d] += myWidth/2;
			else                  myCenter[d] -= myWidth/2;
	}
	Real weight , depth;
	_getSampleDepthAndWeight( densityWeights , temp , position , weightKey , depth , weight );
	depth += depthBias;

	if( depth<minDepth ) depth = Real(minDepth);
	if( depth>maxDepth ) depth = Real(maxDepth);
	int topDepth = int(ceil(depth));

	dx = 1.0-(topDepth-depth);
	if     ( topDepth<=minDepth ) topDepth = minDepth , dx = 1;
	else if( topDepth> maxDepth ) topDepth = maxDepth , dx = 1;

	while( _localDepth( temp )>topDepth ) temp=temp->parent;
	while( _localDepth( temp )<topDepth )
	{
		if( !temp->children ) temp->initChildren( nodeAllocator , _NodeInitializer( *this ) );
		int cIndex = FEMTreeNode::ChildIndex( myCenter , position );
		temp = &temp->children[cIndex];
		myWidth/=2;
		for( int d=0 ; d<Dim ; d++ )
			if( (cIndex>>d) & 1 ) myCenter[d] += myWidth/2;
			else                  myCenter[d] -= myWidth/2;
	}
	width = 1.0 / ( 1<<_localDepth( temp ) );
	_v = v * weight / Real( pow( width , dim ) ) * Real( dx );
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
	_splatPointData< CreateNodes , V >( temp , position , _v , dataInfo , dataKey );
#else // !__GNUC__ || __GNUC__ >=5
	_splatPointData< CreateNodes , V ,  DataSigs ... >( temp , position , _v , dataInfo , dataKey );
#endif // __GNUC__ || __GNUC__ < 4
	if( fabs(1.0-dx) > 1e-6 )
	{
		dx = Real(1.0-dx);
		temp = temp->parent;
		width = 1.0 / ( 1<<_localDepth( temp ) );

		_v = v * weight / Real( pow( width , dim ) ) * Real( dx );
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
		_splatPointData< CreateNodes , V >( temp , position , _v , dataInfo , dataKey );
#else // !__GNUC__ || __GNUC__ >=5
		_splatPointData< CreateNodes , V , DataSigs ... >( temp , position , _v , dataInfo , dataKey );
#endif // __GNUC__ || __GNUC__ < 4
	}
	return weight;
}
template< unsigned int Dim , class Real >
template< bool CreateNodes , unsigned int WeightDegree , class V , unsigned int ... DataSigs >
Real FEMTree< Dim , Real >::_multiSplatPointData( const DensityEstimator< WeightDegree >* densityWeights , FEMTreeNode* node , Point< Real , Dim > position , V v , SparseNodeData< V , UIntPack< DataSigs ... > >& dataInfo , PointSupportKey< IsotropicUIntPack< Dim , WeightDegree > >& weightKey , PointSupportKey< UIntPack< FEMSignature< DataSigs >::Degree ... > >& dataKey , int dim )
{
	typedef UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportSize ... > SupportSizes;
	Real _depth , weight;
	if( densityWeights ) _getSampleDepthAndWeight( *densityWeights , position , weightKey , _depth , weight );
	else weight = (Real)1.;
	V _v = v * weight;

	double values[ Dim ][ SupportSizes::Max() ];
	dataKey.template getNeighbors< CreateNodes >( node , nodeAllocator , _NodeInitializer( *this ) );

	for( FEMTreeNode* _node=node ; _localDepth( _node )>=0 ; _node=_node->parent )
	{
		V __v = _v * (Real)pow( 1<<_localDepth( _node ) , dim );
		Point< Real , Dim > start;
		Real w;
		_startAndWidth( _node , start , w );
		__SetBSplineComponentValues< Real , FEMSignature< DataSigs >::Degree ... >( &position[0] , &start[0] , w , &values[0][0] , SupportSizes::Max() );
		typename FEMTreeNode::template Neighbors< UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportSize ... > >& neighbors = dataKey.neighbors[ _localToGlobal( _localDepth( _node ) ) ];
		double scratch[Dim+1];
		scratch[0] = 1.;
		WindowLoop< Dim >::Run
		(
			ZeroUIntPack< Dim >() , UIntPack< BSplineSupportSizes< FEMSignature< DataSigs >::Degree >::SupportSize ... >() ,
			[&]( int d , int i ){ scratch[d+1] = scratch[d] * values[d][i]; } ,
			[&]( FEMTreeNode* node ){ if( IsActiveNode< Dim >( node ) ) dataInfo[ node ] += __v * (Real)scratch[Dim];	} ,
			neighbors.neighbors()
		);
	}
	return weight;
}

template< unsigned int Dim , class Real >
template< unsigned int WeightDegree , class V , unsigned int ... DataSigs >
Real FEMTree< Dim , Real >::_nearestMultiSplatPointData( const DensityEstimator< WeightDegree >* densityWeights , FEMTreeNode* node , Point< Real , Dim > position , V v , SparseNodeData< V , UIntPack< DataSigs ... > >& dataInfo , PointSupportKey< IsotropicUIntPack< Dim , WeightDegree > >& weightKey , int dim )
{
	Real _depth , weight;
	if( densityWeights ) _getSampleDepthAndWeight( *densityWeights , position , weightKey , _depth , weight );
	else weight = (Real)1.;
	V _v = v * weight;

	for( FEMTreeNode* _node=node ; _localDepth( _node )>=0 ; _node=_node->parent ) if( IsActiveNode< Dim >( _node ) )  dataInfo[ _node ] += _v * (Real)pow( 1<<_localDepth( _node ) , dim );
	return weight;
}
//////////////////////////////////
// MultiThreadedWeightEvaluator //
//////////////////////////////////
template< unsigned int Dim , class Real >
template< unsigned int DensityDegree >
FEMTree< Dim , Real >::MultiThreadedWeightEvaluator< DensityDegree >::MultiThreadedWeightEvaluator( const FEMTree< Dim , Real >* tree , const DensityEstimator< DensityDegree >& density , int threads ) : _density( density ) , _tree( tree )
{
	_threads = std::max< int >( 1 , threads );
	_neighborKeys.resize( _threads );
	for( int t=0 ; t<_neighborKeys.size() ; t++ ) _neighborKeys[t].set( tree->_localToGlobal( density.kernelDepth() ) );
}
template< unsigned int Dim , class Real >
template< unsigned int DensityDegree >
Real FEMTree< Dim , Real >::MultiThreadedWeightEvaluator< DensityDegree >::weight( Point< Real , Dim > p , int thread )
{
	ConstPointSupportKey< IsotropicUIntPack< Dim , DensityDegree > >& nKey = _neighborKeys[thread];
	Real depth , weight;
	_tree->_getSampleDepthAndWeight( _density , p , nKey , depth , weight );
	return weight;
}
