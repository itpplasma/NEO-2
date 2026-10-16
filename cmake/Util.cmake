include(FetchContent)

function(find_or_fetch DEPENDENCY)
    # -D<dep>_SOURCE_DIR=<path> builds against a local checkout; else fetch.
    # FetchContent caches <dep>_SOURCE_DIR as well, so a cached path inside the
    # fetch directory is our own checkout and is refetched when the ref changes.
    string(TOUPPER ${DEPENDENCY} DEPENDENCY_UPPER)
    set(FETCHED_DIR ${FETCHCONTENT_BASE_DIR}/${DEPENDENCY}-src)
    if(${DEPENDENCY}_SOURCE_DIR AND NOT ${DEPENDENCY}_SOURCE_DIR STREQUAL FETCHED_DIR)
        message(STATUS "Using local ${DEPENDENCY} in ${${DEPENDENCY}_SOURCE_DIR}")
    elseif(${DEPENDENCY}_SOURCE_DIR AND (NOT DEFINED ${DEPENDENCY_UPPER}_REF
            OR ${DEPENDENCY_UPPER}_REF STREQUAL ${DEPENDENCY}_FETCHED_REF))
        message(STATUS "Using fetched ${DEPENDENCY} in ${${DEPENDENCY}_SOURCE_DIR}")
    else()
        if(ARGC GREATER 1)
            set(REPO_URL "${ARGV1}")
        else()
            set(REPO_URL https://github.com/itpplasma/${DEPENDENCY}.git)
        endif()

        # Check if a specific ref (branch, tag, or SHA) is provided for this dependency.
        # Pass -DLIBNEO_REF=<ref> on the cmake command line to override.
        if(DEFINED ${DEPENDENCY_UPPER}_REF)
            set(REMOTE_BRANCH ${${DEPENDENCY_UPPER}_REF})
            message(STATUS "Using ${DEPENDENCY} ref ${REMOTE_BRANCH} from ${REPO_URL}")
        else()
            get_branch_or_main(${REPO_URL} REMOTE_BRANCH)
            message(STATUS "Using ${DEPENDENCY} branch ${REMOTE_BRANCH} from ${REPO_URL}")
        endif()

        FetchContent_Declare(
            ${DEPENDENCY}
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
            GIT_REPOSITORY ${REPO_URL}
            GIT_TAG ${REMOTE_BRANCH}
        )
        FetchContent_Populate(${DEPENDENCY})
        set(${DEPENDENCY}_FETCHED_REF ${REMOTE_BRANCH} CACHE INTERNAL
            "${DEPENDENCY} ref of the fetched checkout")
    endif()

    add_subdirectory(${${DEPENDENCY}_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/${DEPENDENCY}
        EXCLUDE_FROM_ALL
    )
endfunction()


function(get_branch_or_main REPO_URL REMOTE_BRANCH)
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    execute_process(
        COMMAND git ls-remote --heads ${REPO_URL} ${BRANCH}
        OUTPUT_VARIABLE BRANCH_EXISTS
        ERROR_VARIABLE GIT_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(BRANCH_EXISTS)
        set(${REMOTE_BRANCH} ${BRANCH} PARENT_SCOPE)
    else()
        message(STATUS "No branch ${BRANCH} at ${REPO_URL}. Return main instead!")
        set(${REMOTE_BRANCH} "main" PARENT_SCOPE)
    endif()
endfunction()
