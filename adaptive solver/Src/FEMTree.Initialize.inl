/*
Copyright (c) 2016, Michael Kazhdan
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

////////////////////////
// FEMTreeInitializer //
////////////////////////
template< unsigned int Dim , class Real >
int FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& node , int maxDepth , std::function< bool ( int , int[] ) > Refine , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
{
	int count = 0;
	int d , off[3];
	node.depthAndOffset( d , off );
	if( node.depth()<maxDepth && Refine( d , off ) )
	{
		node.initChildren( nodeAllocator , NodeInitializer ) , count += 1<<Dim;
		for( int c=0 ; c<(1<<Dim) ; c++ ) count += Initialize( node.children[c] , maxDepth , Refine , nodeAllocator , NodeInitializer );
	}
	return count;
}

template< unsigned int Dim , class Real >
int FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , InputPointStream< Real , Dim >& pointStream , int maxDepth , std::vector< PointSample >& samplePoints , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
{
	auto Leaf = [&]( FEMTreeNode& root , Point< Real , Dim > p , int maxDepth )
	{
		for( int d=0 ; d<Dim ; d++ ) if( p[d]<0 || p[d]>1 ) return (FEMTreeNode*)NULL;
		Point< Real , Dim > center;
		for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
		Real width = Real(1.0);
		FEMTreeNode* node = &root;
		int d = 0;
		while( d<maxDepth )
		{
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
			int cIndex = FEMTreeNode::ChildIndex( center , p );
			node = node->children + cIndex;
			d++;
			width /= 2;
			for( int dd=0 ; dd<Dim ; dd++ )
				if( (cIndex>>dd) & 1 ) center[dd] += width/2;
				else                   center[dd] -= width/2;
		}
		return node;
	};

	// Add the point data
	int outOfBoundPoints = 0 , pointCount = 0;
	{
		std::vector< int > nodeToIndexMap;
		Point< Real , Dim > p;
		while( pointStream.nextPoint( p ) )
		{
			Real weight = (Real)1.;
			FEMTreeNode* temp = Leaf( root , p , maxDepth );
			if( !temp ){ outOfBoundPoints++ ; continue; }
			int nodeIndex = temp->nodeData.nodeIndex;
			if( nodeIndex>=nodeToIndexMap.size() ) nodeToIndexMap.resize( nodeIndex+1 , -1 );
			int idx = nodeToIndexMap[ nodeIndex ];
			if( idx==-1 )
			{
				idx = (int)samplePoints.size();
				nodeToIndexMap[ nodeIndex ] = idx;
				samplePoints.resize( idx+1 ) , samplePoints[idx].node = temp;
			}
			samplePoints[idx].sample += ProjectiveData< Point< Real , Dim > , Real >( p*weight , weight );
			pointCount++;
		}
		pointStream.reset();
	}
	if( outOfBoundPoints  ) WARN( "Found out-of-bound points: %d" , outOfBoundPoints );
	FEMTree< Dim , Real >::MemoryUsage();
	return pointCount;
}

// modified by dojo
// aux 가 없는 conf pts or normalized pts version
template< unsigned int Dim , class Real >
template< class Data >
int FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , InputPointStreamWithData< Real , Dim , Data >& pointStream , int maxDepth , std::vector< PointSample >& samplePoints , std::vector< Data >& sampleData , bool mergeNodeSamples , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer , 
	std::function< Real ( const Point< Real , Dim >&, Data& ) > ProcessData )
{
	auto Leaf = [&]( FEMTreeNode& root , Point< Real , Dim > p , int maxDepth )
	{
		for( int d=0 ; d<Dim ; d++ ) if( p[d]<0 || p[d]>1 ) return (FEMTreeNode*)NULL;
		Point< Real , Dim > center;
		for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
		Real width = Real(1.0);
		FEMTreeNode* node = &root;
		int d = 0;
		while( d<maxDepth )
		{
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
			int cIndex = FEMTreeNode::ChildIndex( center , p );
			node = node->children + cIndex;
			d++;
			width /= 2;
			for( int dd=0 ; dd<Dim ; dd++ )
				if( (cIndex>>dd) & 1 ) center[dd] += width/2;
				else                   center[dd] -= width/2;
		}
		return node;
	};

	// Add the point data
	// 현재 프로세스에서는 이 함수가 ( callback the ProcessData ) 호출된다.
	int outOfBoundPoints = 0 , badData = 0 , pointCount = 0;
	{
		std::vector< int > nodeToIndexMap;
		Point< Real , Dim > p;
		Data d;

		while( pointStream.nextPoint( p , d ) )
		{
			Real weight = ProcessData( p , d );
			if( weight<=0 ){ badData++ ; continue; }
			FEMTreeNode* temp = Leaf( root , p , maxDepth );
			if( !temp ){ outOfBoundPoints++ ; continue; }
			int nodeIndex = temp->nodeData.nodeIndex;
			if( mergeNodeSamples )
			{
				// 같은 node 에 속하는 points 를 merge (summation)
				if( nodeIndex>=nodeToIndexMap.size() ) nodeToIndexMap.resize( nodeIndex+1 , -1 );
				int idx = nodeToIndexMap[ nodeIndex ];
				if( idx==-1 )
				{
					idx = (int)samplePoints.size();
					nodeToIndexMap[ nodeIndex ] = idx;
					samplePoints.resize( idx+1 ) , samplePoints[idx].node = temp;
					sampleData.resize(idx + 1);
				}

				samplePoints[idx].sample += ProjectiveData< Point< Real , Dim > , Real >( p*weight , weight );
				sampleData[ idx ] += d*weight;
			}
			else
			{
				// 모든 points 를 node 에 할당 (여기서는 사용 안 함)
				int idx = (int)samplePoints.size();
				samplePoints.resize( idx+1 ) , sampleData.resize( idx+1 );
				samplePoints[idx].node = temp;
				samplePoints[idx].sample = ProjectiveData< Point< Real , Dim > , Real >( p*weight , weight );
				sampleData[idx] = d * weight;
			}
			pointCount++;
		}
		pointStream.reset();
	}
	if( outOfBoundPoints  ) WARN( "Found out-of-bound points: %d" , outOfBoundPoints );
	if( badData           ) WARN( "Found bad data: %d" , badData );
	FEMTree< Dim , Real >::MemoryUsage();
	return pointCount;
}

// new version (weighted sum 에서 max sample 로...)
// 노이즈 분포 점들의 contribution 을 최소화하기 위함.
// modified by dojo
// aux 가 있는 conf pts version
template< unsigned int Dim, class Real >
template< class Data >
int FEMTreeInitializer< Dim, Real >::Initialize(FEMTreeNode& root, InputPointStreamWithData< Real, Dim, Data >& pointStream, int maxDepth, std::vector< PointSample >& samplePoints, std::vector< Data >& sampleData, bool mergeNodeSamples, Allocator< FEMTreeNode >* nodeAllocator, std::function< void(FEMTreeNode&) > NodeInitializer,
	std::function< Real(const Point< Real, Dim >&, const int, bool&, Data&) > ProcessData)
{
	// modified by dojo
	auto Leaf = [&](FEMTreeNode& root, Point< Real, Dim > p, int maxDepth, bool is_confidence_point)
	{
		for (int d = 0; d < Dim; d++) if (p[d] < 0 || p[d]>1) return (FEMTreeNode*)NULL;
		Point< Real, Dim > center;
		for (int d = 0; d < Dim; d++) center[d] = (Real)0.5;
		Real width = Real(1.0);
		FEMTreeNode* node = &root;
		int d = 0;
		while (d < maxDepth)
		{
			if (!node->children)
			{
				node->initChildren(nodeAllocator, NodeInitializer);
				node->confidence_flag = is_confidence_point ? CONFI_VF_VISITED : CONFI_VF_AUX_ONLY;
			}

			int cIndex = FEMTreeNode::ChildIndex(center, p);
			node = node->children + cIndex;
			if (node->confidence_flag != CONFI_VF_VISITED)
				node->confidence_flag = is_confidence_point ? CONFI_VF_VISITED : CONFI_VF_AUX_ONLY;
			d++;
			width /= 2;
			for (int dd = 0; dd < Dim; dd++)
				if ((cIndex >> dd) & 1) center[dd] += width / 2;
				else                   center[dd] -= width / 2;
		}
		return node;
	};

	std::cout << "Aux. initialization on!" << std::endl;

	// Add the point data
	// modified by dojo
	// 현재 프로세스에서는 이 함수가 ( callback the ProcessData ) 호출된다.
	// mergeNodeSamples == true
	int outOfBoundPoints = 0, badData = 0, pointCount = 0;
	{
		std::vector< int > nodeToIndexMap;
		Point< Real, Dim > p;
		Data d;

		// modified by dojo
		int idx = 0; // note that pointCount is counted only when the point is valid (not bad and inbound)
		std::vector< Real > sample_avr_weights;

		while (pointStream.nextPoint(p, d))
		{
			// modified by dojo
			bool is_confidence_point = true;
			// ProcessData 에서 index 에 기반하여 is_confidence_point 가 재설정
			Real weight = ProcessData(p, idx++, is_confidence_point, d);
			if (weight <= 0) { badData++; continue; }
			// node 의 confidence_flag 설정
			FEMTreeNode* temp = Leaf(root, p, maxDepth, is_confidence_point);

			if (!temp) { outOfBoundPoints++; continue; }
			int nodeIndex = temp->nodeData.nodeIndex;

			// 같은 node 에 속하는 points 를 merge (summation)
			// node 를 cell 로 생각하면 됨.
			if (nodeIndex >= nodeToIndexMap.size()) nodeToIndexMap.resize(nodeIndex + 1, -1);
			int idx = nodeToIndexMap[nodeIndex];
			bool is_new = false;
			if (idx == -1)
			{
				idx = (int)samplePoints.size();
				nodeToIndexMap[nodeIndex] = idx;
				samplePoints.resize(idx + 1), samplePoints[idx].node = temp;
				sampleData.resize(idx + 1);
				// modified by dojo
				sample_avr_weights.resize(idx + 1), sample_avr_weights[idx] = 0;
				//is_new = true;
			}

			//if (is_new)
			//{
			//	samplePoints[idx].sample = ProjectiveData< Point< Real, Dim >, Real >(p*weight, weight);
			//	sampleData[idx] = d * weight;
			//}
			//else
			{
				// temp (leaf) node 에 conf 가 한 번이라도 visit 했으면 'CONFI_VF_VISITED'
				// 여기에는 aux point 가 영향 안 주도록 한다. (&& is_confidence_point)
				//if (temp->confidence_flag != CONFI_VF_AUX_ONLY) // which means temp is CONFI_VF_VISITED
				//{
					//if (is_confidence_point)
					{
					// weight is 1
					// 그리고 aux 가 conf pts 다음에 들어 오므로, sample 을 reset 할 필요 없다.
					samplePoints[idx].sample += ProjectiveData< Point< Real, Dim >, Real >(p*weight, weight);
					sampleData[idx] += d * weight;
					sample_avr_weights[idx] += weight * weight;
					}
				//}
				//else
				//{
				//	// temp->confidence_flag == CONFI_VF_AUX_ONLY 이면 !is_confidence_point 이다
				//	ProjectiveData< Point< Real, Dim >, Real >& sample = samplePoints[idx].sample;
				//	if (sample.weight < weight)
				//	{
				//		samplePoints[idx].sample = ProjectiveData< Point< Real, Dim >, Real >(p*weight, weight);
				//		sampleData[idx] = d * weight;
				//		sample_avr_weights[idx] = weight * weight;
				//	}
				//	//samplePoints[idx].sample += ProjectiveData< Point< Real, Dim >, Real >(p*weight, weight);
				//	//sampleData[idx] += d * weight;
				//	//sample_avr_weights[idx] += weight * weight;
				//}
			}
			pointCount++;
		}
		pointStream.reset();

		// modified by dojo //
#pragma omp parallel for
		for (int i = 0; i < samplePoints.size(); i++)
		{
			//if (samplePoints[i].node->confidence_flag == CONFI_VF_AUX_ONLY)
			{
				ProjectiveData< Point< Real, Dim >, Real >& sample = samplePoints[i].sample;
				Real w_sum = sample.weight;
				Real w_avr = sample_avr_weights[i] / w_sum;
				sampleData[i] /= w_sum;

				Point< Real, Dim > p_avr = sample.data / w_sum;
				samplePoints[i].sample = ProjectiveData< Point< Real, Dim >, Real >(p_avr * w_avr, w_avr);
				//if (samplePoints[i].node->confidence_flag == CONFI_VF_VISITED
				//	&& w_avr < 0.99) printf("GGG << %f\n", w_avr);
				sampleData[i] *= w_avr;
			}
		}
	}
	if (outOfBoundPoints) WARN("Found out-of-bound points: %d", outOfBoundPoints);
	if (badData) WARN("Found bad data: %d", badData);

	FEMTree< Dim, Real >::MemoryUsage();
	return pointCount;
}

template< unsigned int Dim , class Real >
void FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , const std::vector< Point< Real , Dim > >& vertices , const std::vector< SimplexIndex< Dim-1 > >& simplices , int maxDepth , std::vector< PointSample >& samples , bool mergeNodeSamples , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
{
	std::vector< int > nodeToIndexMap;
#pragma omp parallel for
	for( int i=0 ; i<simplices.size() ; i++ )
	{
		Simplex< Real , Dim , Dim-1 > s;
		for( int k=0 ; k<Dim ; k++ ) s[k] = vertices[ simplices[i][k] ];
		int sCount;
		if( mergeNodeSamples ) sCount = _AddSimplex( root , s , maxDepth , samples , &nodeToIndexMap , nodeAllocator , NodeInitializer );
		else                   sCount = _AddSimplex( root , s , maxDepth , samples , NULL ,            nodeAllocator , NodeInitializer );
	}
	FEMTree< Dim , Real >::MemoryUsage();
}

template< unsigned int Dim , class Real >
int FEMTreeInitializer< Dim , Real >::_AddSimplex( FEMTreeNode& root , Simplex< Real , Dim , Dim-1 >& s , int maxDepth , std::vector< PointSample >& samples , std::vector< int >* nodeToIndexMap , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
{
	std::vector< Simplex< Real , Dim , Dim-1 > > subSimplices;
	subSimplices.push_back( s );

	// Clip the simplex to the unit cube
	{
		for( int d=0 ; d<Dim ; d++ )
		{
			Point< Real , Dim > n;
			n[d] = 1;
			{
				std::vector< Simplex< Real , Dim , Dim-1 > > back , front;
				for( int i=0 ; i<subSimplices.size() ; i++ ) subSimplices[i].split( n , 0 , back , front );
				subSimplices = front;
			}
			{
				std::vector< Simplex< Real , Dim , Dim-1 > > back , front;
				for( int i=0 ; i<subSimplices.size() ; i++ ) subSimplices[i].split( n , 1 , back , front );
				subSimplices = back;
			}
		}
	}

	struct RegularGridIndex
	{
		int idx[Dim];
		bool operator != ( const RegularGridIndex& i ) const
		{
			for( int d=0 ; d<Dim ; d++ ) if( idx[d]!=i.idx[d] ) return true;
			return false;
		}
	};

	auto Leaf = [&]( Point< Real , Dim > p , int maxDepth )
	{
		for( int d=0 ; d<Dim ; d++ ) if( p[d]<0 || p[d]>1 ) return (FEMTreeNode*)NULL;
		Point< Real , Dim > center;
		for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
		Real width = Real(1.0);
		FEMTreeNode* node = &root;
		int d=0;
		while( d<maxDepth )
		{
#pragma omp critical
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
			int cIndex = FEMTreeNode::ChildIndex( center , p );
			node = node->children + cIndex;
			d++;
			width /= 2;
			for( int d=0 ; d<Dim ; d++ )
				if( (cIndex>>d) & 1 ) center[d] += width/2;
				else                  center[d] -= width/2;
		}
		return node;
	};


	int sCount = 0;
	for( int i=0 ; i<subSimplices.size() ; i++ )
	{
		// Find the finest depth at which the simplex is entirely within a node
		int tDepth;
		RegularGridIndex idx0 , idx;
		for( tDepth=0 ; tDepth<maxDepth ; tDepth++ )
		{
			// Get the grid index of the first vertex of the simplex
			for( int d=0 ; d<Dim ; d++ ) idx0.idx[d] = idx.idx[d] = (int)( subSimplices[i][0][d] * (1<<(tDepth+1)) );
			bool done = false;
			for( int k=0 ; k<Dim && !done ; k++ )
			{
				for( int d=0 ; d<Dim ; d++ ) idx.idx[d] = (int)( subSimplices[i][k][d] * (1<<(tDepth+1)) );
				if( idx!=idx0 ) done = true;
			}
			if( done ) break;
		}

		// Generate a point in the middle of the simplex
		for( int i=0 ; i<subSimplices.size() ; i++ ) sCount += _AddSimplex( Leaf( subSimplices[i].center() , tDepth ) , subSimplices[i] , maxDepth , samples , nodeToIndexMap , nodeAllocator , NodeInitializer );
	}
	return sCount;
}
template< unsigned int Dim , class Real >
int FEMTreeInitializer< Dim , Real >::_AddSimplex( FEMTreeNode* node , Simplex< Real , Dim , Dim-1 >& s , int maxDepth , std::vector< PointSample >& samples , std::vector< int >* nodeToIndexMap , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
{
	int d = node->depth();
	if( d==maxDepth )
	{
		Real weight = s.measure();
		Point< Real , Dim > position = s.center() , normal;
		{
			Point< Real , Dim > v[Dim-1];
			for( int k=0 ; k<Dim-1 ; k++ ) v[k] = s[k+1]-s[0];
			normal = Point< Real , Dim >::CrossProduct( v );
		}
		if( weight && weight==weight )
		{
			if( nodeToIndexMap )
			{
				int nodeIndex = node->nodeData.nodeIndex;
#pragma omp critical
				{
					if( nodeIndex>=nodeToIndexMap->size() ) nodeToIndexMap->resize( nodeIndex+1 , -1 );
					int idx = (*nodeToIndexMap)[ nodeIndex ];
					if( idx==-1 )
					{
						idx = (int)samples.size();
						(*nodeToIndexMap)[ nodeIndex ] = idx;
						samples.resize( idx+1 );
						samples[idx].node = node;
					}
					samples[idx].sample += ProjectiveData< Point< Real , Dim > , Real >( position*weight , weight );
				}
			}
			else
			{
#pragma omp critical
				{
					int idx = (int)samples.size();
					samples.resize( idx+1 );
					samples[idx].node = node;
					samples[idx].sample = ProjectiveData< Point< Real , Dim > , Real >( position*weight , weight );
				}
			}
		}
		return 1;
	}
	else
	{
		int sCount = 0;
#pragma omp critical
		if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );

		// Split up the simplex and pass the parts on to the children
		Point< Real , Dim > center;
		Real width;
		node->centerAndWidth( center , width );

		std::vector< std::vector< Simplex< Real , Dim , Dim-1 > > > childSimplices( 1 );
		childSimplices[0].push_back( s );
		for( int d=0 ; d<Dim ; d++ )
		{
			Point< Real , Dim > n ; n[Dim-d-1] = 1;
			std::vector< std::vector< Simplex< Real , Dim , Dim-1 > > > temp( (int)( 1<<(d+1) ) );
			for( int c=0 ; c<(1<<d) ; c++ ) for( int i=0 ; i<childSimplices[c].size() ; i++ ) childSimplices[c][i].split( n , center[Dim-d-1] , temp[2*c] , temp[2*c+1] );
			childSimplices = temp;
		}
		for( int c=0 ; c<(1<<Dim) ; c++ ) for( int i=0 ; i<childSimplices[c].size() ; i++ ) if( childSimplices[c][i].measure() ) sCount += _AddSimplex( node->children+c , childSimplices[c][i] , maxDepth , samples , nodeToIndexMap , nodeAllocator , NodeInitializer );
		return sCount;
	}
}

template< unsigned int Dim , class Real >
void FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , const std::vector< Point< Real , Dim > >& vertices , const std::vector< SimplexIndex< Dim-1 > >& simplices , int maxDepth , std::vector< NodeSimplices< Dim , Real > >& nodeSimplices , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
{
	std::vector< int > nodeToIndexMap;
	for( int i=0 ; i<simplices.size() ; i++ )
	{
		Simplex< Real , Dim , Dim-1 > s;
		for( int k=0 ; k<Dim ; k++ ) s[k] = vertices[ simplices[i][k] ];
		_AddSimplex( root , s , maxDepth , nodeSimplices , nodeToIndexMap , nodeAllocator , NodeInitializer );
	}
	FEMTree< Dim , Real >::MemoryUsage();
}

template< unsigned int Dim , class Real >
int FEMTreeInitializer< Dim , Real >::_AddSimplex( FEMTreeNode& root , Simplex< Real , Dim , Dim-1 >& s , int maxDepth , std::vector< NodeSimplices< Dim , Real > >& simplices , std::vector< int >& nodeToIndexMap , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
{
	std::vector< Simplex< Real , Dim , Dim-1 > > subSimplices;
	subSimplices.push_back( s );

	// Clip the simplex to the unit cube
	{
		for( int d=0 ; d<Dim ; d++ )
		{
			Point< Real , Dim > n;
			n[d] = 1;
			{
				std::vector< Simplex< Real , Dim , Dim-1 > > back , front;
				for( int i=0 ; i<subSimplices.size() ; i++ ) subSimplices[i].split( n , 0 , back , front );
				subSimplices = front;
			}
			{
				std::vector< Simplex< Real , Dim , Dim-1 > > back , front;
				for( int i=0 ; i<subSimplices.size() ; i++ ) subSimplices[i].split( n , 1 , back , front );
				subSimplices = back;
			}
		}
	}

	struct RegularGridIndex
	{
		int idx[Dim];
		bool operator != ( const RegularGridIndex& i ) const
		{
			for( int d=0 ; d<Dim ; d++ ) if( idx[d]!=i.idx[d] ) return true;
			return false;
		}
	};

	auto Leaf = [&]( Point< Real , Dim > p , int maxDepth )
	{
		for( int d=0 ; d<Dim ; d++ ) if( p[d]<0 || p[d]>1 ) return (FEMTreeNode*)NULL;
		Point< Real , Dim > center;
		for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
		Real width = Real(1.0);
		FEMTreeNode* node = &root;
		int d=0;
		while( d<maxDepth )
		{
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
			int cIndex = FEMTreeNode::ChildIndex( center , p );
			node = node->children + cIndex;
			d++;
			width /= 2;
			for( int d=0 ; d<Dim ; d++ )
				if( (cIndex>>d) & 1 ) center[d] += width/2;
				else                  center[d] -= width/2;
		}
		return node;
	};

	int sCount = 0;

	for( int i=0 ; i<subSimplices.size() ; i++ )
	{
		// Find the finest depth at which the simplex is entirely within a node
		int tDepth;
		RegularGridIndex idx0 , idx;
		for( tDepth=0 ; tDepth<maxDepth ; tDepth++ )
		{
			// Get the grid index of the first vertex of the simplex
			for( int d=0 ; d<Dim ; d++ ) idx0.idx[d] = (int)( subSimplices[i][0][d] * (1<<(tDepth+1)) );
			bool done = false;
			for( int k=0 ; k<Dim && !done ; k++ )
			{
				for( int d=0 ; d<Dim ; d++ ) idx.idx[d] = (int)( subSimplices[i][k][d] * (1<<(tDepth+1)) );
				if( idx!=idx0 ) done = true;
			}
			if( done ) break;
		}

		// Add the simplex to the node
		FEMTreeNode* subSimplexNode = Leaf( subSimplices[i].center() , tDepth );
		for( int i=0 ; i<subSimplices.size() ; i++ ) sCount += _AddSimplex( subSimplexNode , subSimplices[i] , maxDepth , simplices , nodeToIndexMap , nodeAllocator , NodeInitializer );
	}
	return sCount;
}
template< unsigned int Dim , class Real >
int FEMTreeInitializer< Dim , Real >::_AddSimplex( FEMTreeNode* node , Simplex< Real , Dim , Dim-1 >& s , int maxDepth , std::vector< NodeSimplices< Dim , Real > >& simplices , std::vector< int >& nodeToIndexMap , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
{
	int d = node->depth();
	if( d==maxDepth )
	{
		// If the simplex has non-zero size, add it to the list
		Real weight = s.measure();
		if( weight && weight==weight )
		{
			int nodeIndex = node->nodeData.nodeIndex;
			if( nodeIndex>=nodeToIndexMap.size() ) nodeToIndexMap.resize( nodeIndex+1 , -1 );
			int idx = nodeToIndexMap[ nodeIndex ];
			if( idx==-1 )
			{
				idx = (int)simplices.size();
				nodeToIndexMap[ nodeIndex ] = idx;
				simplices.resize( idx+1 );
				simplices[idx].node = node;
			}
			simplices[idx].data.push_back( s );
		}
		return 1;
	}
	else
	{
		int sCount = 0;
		if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );

		// Split up the simplex and pass the parts on to the children
		Point< Real , Dim > center;
		Real width;
		node->centerAndWidth( center , width );

		std::vector< std::vector< Simplex< Real , Dim , Dim-1 > > > childSimplices( 1 );
		childSimplices[0].push_back( s );
		for( int d=0 ; d<Dim ; d++ )
		{
			Point< Real , Dim > n ; n[Dim-d-1] = 1;
			std::vector< std::vector< Simplex< Real , Dim , Dim-1 > > > temp( (int)( 1<<(d+1) ) );
			for( int c=0 ; c<(1<<d) ; c++ ) for( int i=0 ; i<childSimplices[c].size() ; i++ ) childSimplices[c][i].split( n , center[Dim-d-1] , temp[2*c] , temp[2*c+1] );
			childSimplices = temp;
		}
		for( int c=0 ; c<(1<<Dim) ; c++ ) for( int i=0 ; i<childSimplices[c].size() ; i++ ) sCount += _AddSimplex( node->children+c , childSimplices[c][i] , maxDepth , simplices , nodeToIndexMap , nodeAllocator , NodeInitializer );
		return sCount;
	}
}

template< unsigned int Dim , class Real >
template< class Data , class _Data , bool Dual >
int FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , ConstPointer( Data ) values , ConstPointer( int ) labels , int resolution[Dim] , std::vector< NodeSample< Dim , _Data > > derivatives[Dim] , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer , std::function< _Data ( const Data& ) > DataConverter )
{
	auto Leaf = [&]( FEMTreeNode& root , const int idx[Dim] , int maxDepth )
	{
		for( int d=0 ; d<Dim ; d++ ) if( idx[d]<0 || idx[d]>=(1<<maxDepth) ) return (FEMTreeNode*)NULL;
		FEMTreeNode* node = &root;
		for( int d=0 ; d<maxDepth ; d++ )
		{
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
			int cIndex = 0;
			for( int dd=0 ; dd<Dim ; dd++ ) if( idx[dd]&(1<<(maxDepth-d-1)) ) cIndex |= 1<<dd;
			node = node->children + cIndex;
		}
		return node;
	};
	auto FactorIndex = []( size_t i , const int resolution[Dim] , int idx[Dim] )
	{
		size_t ii = i;
		for( int d=0 ; d<Dim ; d++ ) idx[d] = ii % resolution[d] , ii /= resolution[d];
	};
	auto MakeIndex = [] ( const int idx[Dim] , const int resolution[Dim] )
	{
		size_t i = 0;
		for( int d=0 ; d<Dim ; d++ ) i = i * resolution[Dim-1-d] + idx[Dim-1-d];
		return i;
	};


	int maxResolution = resolution[0];
	for( int d=1 ; d<Dim ; d++ ) maxResolution = std::max< int >( maxResolution , resolution[d] );
	int maxDepth = 0;
	while( ( (1<<maxDepth) + ( Dual ? 0 : 1 ) )<maxResolution ) maxDepth++;

	size_t totalRes = 1;
	for( int d=0 ; d<Dim ; d++ ) totalRes *= resolution[d];

	// Iterate over each direction
	for( int d=0 ; d<Dim ; d++ ) for( size_t i=0 ; i<totalRes ; i++ )
	{
		// Factor the index into directional components and get the index of the next cell
		int idx[Dim] ; FactorIndex( i , resolution , idx ) ; idx[d]++;

		if( idx[d]<resolution[d] )
		{
			// Get the index of the next cell
			size_t ii = MakeIndex( idx , resolution );

			// [NOTE] There are no derivatives across negative labels
			if( labels[i]!=labels[ii] && labels[i]>=0 && labels[ii]>=0 )
			{
				if( !Dual ) idx[d]--;
				NodeSample< Dim , _Data > nodeSample;
				nodeSample.node = Leaf( root , idx , maxDepth );
				nodeSample.data = DataConverter( values[ii] ) - DataConverter( values[i] );
				if( nodeSample.node ) derivatives[d].push_back( nodeSample );
			}
		}
	}
	return maxDepth;
}

template< unsigned int Dim , class Real >
template< bool Dual , class Data >
unsigned int FEMTreeInitializer< Dim , Real >::Initialize( FEMTreeNode& root , DerivativeStream< Data >& dStream , std::vector< NodeSample< Dim , Data > > derivatives[Dim] , Allocator< FEMTreeNode >* nodeAllocator , std::function< void ( FEMTreeNode& ) > NodeInitializer )
{
	// Note:
	// --   Dual: The difference between [i] and [i+1] is stored at cell [i+1]
	// -- Primal: The difference between [i] and [i+1] is stored at cell [i]

	// Find the leaf containing the specified cell index
	auto Leaf = [&]( FEMTreeNode& root , const unsigned int idx[Dim] , unsigned int maxDepth )
	{
		for( int d=0 ; d<Dim ; d++ ) if( idx[d]<0 || idx[d]>=(unsigned int)(1<<maxDepth) ) return (FEMTreeNode*)NULL;
		FEMTreeNode* node = &root;
		for( unsigned int d=0 ; d<maxDepth ; d++ )
		{
			if( !node->children ) node->initChildren( nodeAllocator , NodeInitializer );
			int cIndex = 0;
			for( int dd=0 ; dd<Dim ; dd++ ) if( idx[dd]&(1<<(maxDepth-d-1)) ) cIndex |= 1<<dd;
			node = node->children + cIndex;
		}
		return node;
	};

	unsigned int resolution[Dim];
	dStream.resolution( resolution );
	unsigned int maxResolution = resolution[0];
	for( int d=1 ; d<Dim ; d++ ) maxResolution = std::max< unsigned int >( maxResolution , resolution[d] );
	unsigned int maxDepth = 0;

	// If we are using a dual formulation, we need at least maxResolution cells.
	// Otherwise, we need at least maxResolution-1 cells.
	while( (unsigned int)( (1<<maxDepth) + ( Dual ? 0 : 1 ) )<maxResolution ) maxDepth++;

	unsigned int idx[Dim] , dir;
	Data dValue;
	while( dStream.nextDerivative( idx , dir , dValue ) )
	{
		if( Dual ) idx[dir]++;
		NodeSample< Dim , Data > nodeSample;
		nodeSample.node = Leaf( root , idx , maxDepth );
		nodeSample.data = dValue;
		if( nodeSample.node ) derivatives[dir].push_back( nodeSample );
	}
	return maxDepth;
}
