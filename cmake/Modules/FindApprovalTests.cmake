include(FetchContent)

FetchContent_Declare(
        ApprovalTests
        GIT_REPOSITORY https://github.com/approvals/ApprovalTests.cpp.git
        GIT_TAG v.10.13.0
)
FetchContent_MakeAvailable(ApprovalTests)
