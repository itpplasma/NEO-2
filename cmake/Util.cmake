include(FetchContent)

function(find_or_fetch DEPENDENCY)
    # -D<dep>_SOURCE_DIR=<path> builds against a local checkout; else fetch.
    string(TOUPPER ${DEPENDENCY} DEPENDENCY_UPPER)
    if(${DEPENDENCY}_SOURCE_DIR)
        message(STATUS "Using local ${DEPENDENCY} in ${${DEPENDENCY}_SOURCE_DIR}")
    else()
        set(REPO_URL https://github.com/itpplasma/${DEPENDENCY}.git)

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
