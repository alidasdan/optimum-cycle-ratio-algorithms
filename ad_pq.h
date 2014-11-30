//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
#ifndef AD_PQ_INCLUDED
#define AD_PQ_INCLUDED

#include "ad_globals.h"

// A min priority queue implementation based on binary heaps.

template< class key_t, class info_t >
class ad_pq {
public:

    // Structure of a heap node.
    struct ad_pq_node {
        key_t  key;   // Ordering is on this key.
        info_t info;  // Extra info.
        int pos;      // If heap[i] points to this node, pos is i.

        // Constructors:
        ad_pq_node()
        {        
            /*EMPTY*/
        }
        ad_pq_node( const ad_pq_node& rhs )
        {
            key = rhs.key;
            info = rhs.info;
            pos = rhs.pos;
        }
        ad_pq_node( const key_t& k, const info_t& i, int p )
        {
            key = k;
            info = i;
            pos = p;
        }

        ad_pq_node& operator=( const ad_pq_node& rhs )
        {
            if ( this != &rhs ) {
                key = rhs.key;
                info = rhs.info;
                pos = rhs.pos;
            }
            return *this;
        }
    };

    typedef ad_pq_node* ad_pq_node_ptr;

    // Constructors.
    ad_pq()
    {
        hsize = 0;
        len = MAX_ALLOC_SIZE;
        heap = new ad_pq_node_ptr[ len ];
    }

    ad_pq( int n )
    {
        hsize = 0;
        len = n;
        heap = new ad_pq_node_ptr[ len ];
    }

    // Destructor.
    ~ad_pq()
    {
        for ( int i = 1; i <= hsize; ++i )
            delete heap[ i ];
        delete [] heap;
    }

    // Insert (key, info) into the heap.
    ad_pq_node_ptr put( const key_t& key, const info_t& info )
    {
        ++hsize;

        // Resize heap if it does not have enough space.
        if ( len <= hsize )
            resize();

        // Move nodes down until a suitable place for (key, info) is
        // found. This loop is similar to sift_up but the latter uses
        // exchanges rather than copies.
        int i = hsize;
        while ( ( i > 1 ) && ( key < heap[ parent( i ) ]->key ) ) {
            int p = parent( i );
            copy( i, p );
            i = p;
        }
        heap[ i ] = new ad_pq_node( key, info, i );

        return heap[ i ];
    }  // put

    // Update node ptr's key.
    void update_key( const key_t& key, ad_pq_node_ptr ptr )
    {
        if ( key < ptr->key ) {
            ptr->key = key;
            sift_up( ptr->pos );
        } else if ( key > ptr->key ) {
            ptr->key = key;
            sift_down( ptr->pos );
        }
        /* else do nothing */
    }  // update_key

    // Update node ptr's all information.
    void update_node( const key_t& key, const info_t& info, ad_pq_node_ptr ptr )
    {
        ptr->info = info;
        update_key( key, ptr );
    }

    // Return the properties of the min node:
    key_t getkey() 
    {
#ifdef DEBUG
        assert( hsize > 0 );
#endif
        return heap[ 1 ]->key;
    }
    info_t getinfo()
    {
#ifdef DEBUG
        assert( hsize > 0 );
#endif
        return heap[ 1 ]->info;
    }

#if 0
    // The following functions are not used but they work correctly:

    // Delete the node with the min key.
    void del_min()
    {
        delete heap[ 1 ];
        copy( 1, hsize-- );
        sift_down( 1 );
    }

    // Delete the node ptr.
    void del( ad_pq_node_ptr ptr )
    {
        int i = ptr->pos;
        delete heap[ i ];
        copy( i, hsize-- );
        sift_down( i );
    }

    bool empty()
    {
        return ( hsize == 0 );
    }

    int size()
    {
        return hsize;
    }
#endif

private:
    int len;         // Size of array holding heap.
    int hsize;       // Number of nodes in heap.

    // Nodes are pointed to by pointers in heap[1], heap[2], ...,
    // heap[hsize], where heap[1] points to the root of the heap, the
    // min node in the heap. Note that len must be larger than
    // hsize. Also note that the heap property implies heap[i] <=
    // min(heap[left(i)], heap[right(i)]), or heap[parent(i)] <=
    // heap[i].
    ad_pq_node_ptr *heap;

private:
    // Check if i is valid.
    void check_inx( int i ) const
    {
        assert( ( 1 <= i ) && ( i <= hsize ) );
    }

    // Return parent, left child, or right child of node i:
    int parent( int i ) const
    {
        return ( i >> 1 );
    }
    int left( int i ) const
    {
        return ( i << 1 );
    }
    int right( int i ) const
    {
        return ( i << 1 ) + 1;
    }

    // Copy node j to node i.
    void copy( int i, int j ) 
    {
        heap[ i ] = heap[ j ];
        heap[ i ]->pos = i;
    }

    // Exchange node j with node i.
    void exchange( int i, int j )
    {
        heap[ 0 ] = heap[ i ];
        copy( i, j );
        copy( j, 0 );
    }

    // Move up thru exchanges starting from node i.
    void sift_up( int i )
    {
#ifdef DEBUG
        check_inx( i );
#endif
        int p = parent( i );
        while ( ( i > 1 ) && ( heap[ i ]->key < heap[ p ]->key ) ) {
            exchange( i, p );
            i = p;
            p = parent( p ); // i can also be used as an argument.
        }
    }  // sift_up

    // Move down thru exchanges starting from node i.
    void sift_down( int i )
    {
#ifdef DEBUG
        check_inx( i );
#endif
        int l;
        while ( l = left( i ), l <= hsize ) {
            // smallest = min(heap[i], heap[left(i)], heap[right(i)]).
            int smallest = ( ( heap[ l ]->key < heap[ i ]->key ) ? l : i );
            int r = right( i );
            if ( ( r <= hsize ) && ( heap[ r ]->key < heap[ smallest ]->key ) )
                smallest = r;

            if ( smallest != i ) {
                exchange( i, smallest );
                i = smallest;
            }
            else
                break;
        }
    }  // sift_down

    // Resize heap's size by MAX_ALLOC_SIZE.
    void resize()
    {
        len += MAX_ALLOC_SIZE;
        ad_pq_node_ptr *new_heap = new ad_pq_node_ptr[ len ];
        for ( int j = 1; j < hsize; ++j )
            new_heap[ j ] = heap[ j ];
        delete [] heap;
        heap = new_heap;
    }  // resize
};  // ad_pq

#endif
