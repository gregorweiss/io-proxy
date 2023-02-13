# cmake/ProcessOptions.cmake

# Here all user defined options will be processed.

function(PROCESS_WITH_SIONLIB)
    set(HAVE_SIONLIB OFF)
    if (with-sionlib)
        if (NOT ${with-sionlib} STREQUAL "ON")
            set(SIONLIB_ROOT_DIR "${with-sionlib}" CACHE INTERNAL "sionlib")
        endif ()

        find_package(SIONlib)
        include_directories(${SIONLIB_INCLUDE})

        if (SIONLIB_FOUND)
            set(HAVE_SIONLIB ON CACHE INTERNAL "sionlib")
            add_definitions( -DHAVE_SIONLIB )
        endif ()
    endif ()
endfunction()
