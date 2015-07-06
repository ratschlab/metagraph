#ifndef __UNIX_TOOLS_HPP__
#define __UNIX_TOOLS_HPP__

/**
 * This header file collects some useful functions to find out about
 * system resources and stuff
 */

#include <stdio.h>

/** This returns the currently used memory by the process.
 *
 * The code was copied and has been modified from:
 * http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use
 */
size_t
get_curr_mem () {
    long rss = 0L;
    FILE* fp = NULL;

    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
        return (size_t)0L;      /* Can't open? */

    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
        fclose( fp );
        return (size_t)0L;      /* Can't read? */
    }

    fclose( fp );
    return (size_t)rss * (size_t) sysconf( _SC_PAGESIZE);
}


#endif
