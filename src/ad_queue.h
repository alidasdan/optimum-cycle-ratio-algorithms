//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
#ifndef AD_QUEUE_INCLUDED
#define AD_QUEUE_INCLUDED

// Queue implementation.

template< class T >
class ad_queue {
public:

    // Constructors:
    ad_queue()
    {
        init();
        len = 0;
    }

    ad_queue( int l )
    {
        init();
        len = l;
        buf = new T[ len ];
    }

    // Destructor.
    ~ad_queue()
    {
        delete [] buf;
    }

    int length() const
    {
        return len;
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
    }  // put

    T get()
    {
#ifdef DEBUG
        assert( is_not_empty() );
        cur_len--;
#endif
        return buf[ head++ ];
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
        return ( tail >= len );
    }
    bool is_not_full() const
    {
        // = ! is_full().
        return ( tail < len );
    }

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
