from libc.stdlib cimport malloc, free
from libc.string cimport strdup, strcpy

cdef extern from "find_pair.c":
    int find_pair_entry(int argc, char **argv)

def run_find_pair(list args):
    """
    Runs the find_pair CLI equivalent from Python
    Example:
        run_find_pair([])
    """
    cdef int argc = <Py_ssize_t>(len(args))
    cdef char **argv = <char **>malloc(<size_t>(argc * sizeof(char *)))
    cdef int i, status
    cdef int length
    cdef bytes s
    cdef char *c_str

    if not argv:
        raise MemoryError("Cannot allocate array")

    for i in range(argc):
        s = args[i].encode("utf-8")  # get bytes representation
        length = len(s)
        c_str = <char *>malloc(<size_t>(length + 1))  # +1 for null terminator
        if not c_str:
            # Free already allocated strings before raising
            for j in range(i):
                free(<void *>argv[j])
            free(<void *>argv)
            raise MemoryError("Cannot allocate string buffer")
        strcpy(c_str, s)    # copy string contents (C expects null-terminated)
        (<char **>argv)[i] = c_str

    status = find_pair_entry(argc, argv)

    for i in range(argc):
        if argv[i] != NULL:
            free(<void *>argv[i])
    free(<void *>argv)

    return status