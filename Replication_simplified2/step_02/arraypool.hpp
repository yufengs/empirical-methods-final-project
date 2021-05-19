#pragma once


class ArrayPool {
    char * base;
    size_t size;
    size_t used;
public:
    ArrayPool()
    {
        base = NULL;
        size = 0;
        used = 0;
    }

    ~ArrayPool()
    {
        mexPrintf("ArrayPool::destructor\n");
        clear();
    }

    void clear()
    {
        free(base);
        base = NULL;
        size = 0;
        used = 0;
    }

    void resize(size_t n)
    {
        if(n != size)
        {
            mexPrintf("ArrayPool::alloc %d\n", n);

            clear();
            base = (char*)malloc(n);
            size = n;
            used = 0;
        }
        else
        {
            used = 0;
        }
    }

    void * alloc(size_t n)
    {
        if( (used+n) > size )
        {
            mexPrintf("excedded array pool size\n");
            return NULL;
        }

        void * ptr = (void *) ((uintptr_t)base + used);

        used += n;

        return ptr;
    }
};

#if _WIN32
extern __declspec(thread) ArrayPool * arrayPool;
#else
extern __thread ArrayPool * arrayPool;
#endif
