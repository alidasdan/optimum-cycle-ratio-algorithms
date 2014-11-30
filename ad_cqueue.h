//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
#ifndef AD_CQUEUE_INCLUDED
#define AD_CQUEUE_INCLUDED

// A circular queue implementation.

template< class T >
class ad_cqueue {
public:

    // Constructors:
    ad_cqueue()
    {
        init();
        len = 1;
    }

    ad_cqueue( int l )
    {
        init();
        len = l + 1;
        buf = new T[ len ];
    }

    // Destructor.
    ~ad_cqueue()
    {
        delete [] buf;
    }

    int length() const
    {
        // The interface does not reveal the extra location.
        return len - 1;
    }

#ifdef DEBUG
    int level()
    {
        return cur_len;
    }
#endif      

    void init()
    {
        head = tail = 0;
#ifdef DEBUG
        cur_len = 0;
#endif
    }

    void init_head()
    {
        head = 0;
    }

    void put( int v )
    {
#ifdef DEBUG
        assert( is_not_full() );
        cur_len++;
#endif
        buf[ tail++ ] = v;

        // tail = ( tail + 1 ) % len;
        if ( tail == len )
            tail = 0;
    }  // put

    T get()
    {
#ifdef DEBUG
        assert( is_not_empty() );
        cur_len--;
#endif
        T v = buf[ head++ ];

        // head = ( head + 1 ) % len;
        if ( head == len )
            head = 0;
        return v;
    }  // get

    bool is_empty() const
    {
        return ( head == tail );
    }
    bool is_not_empty() const
    {
        // = ! is_empty().
        return ( head != tail );
    }
    bool is_full() const
    {
        return ( head == ( tail + 1 ) % len );
    }
    bool is_not_full() const
    {
        // = ! is_full().
        return ( head != ( tail + 1 ) % len );
    }

    void print() const
    {
        int j = 0;
        printf( "printing queue elements - start:\n" );
        for ( int i = head; i != tail; i = ( i + 1 ) % len ) {
            printf( "queue[ %d ] = %d\n", i, buf[ i ] + 1 );
            j++;
            assert( j <= len );
        }
        printf( "printing queue elements - end:\n" );
    }  // print

private:
    int len;
    int head;
    int tail;

#ifdef DEBUG
    int cur_len;  // Current number of elements, queried using level().
#endif

    // Elements are in buf[ head ], buf[ head+1 ], ..., buf[ tail - 1 ].
    // The next location to put an element into is buf[ tail ].
    T *buf;
};

#endif
