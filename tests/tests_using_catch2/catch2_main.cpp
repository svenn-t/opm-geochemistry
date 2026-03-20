// This tells Catch to provide a main() - only do this in one cpp file!
#if __has_include(<catch2/catch.hpp>)
#define APPROVALS_CATCH
#else
#define APPROVALS_CATCH2_V3
#endif
#include "ApprovalTests.hpp"

// This puts "received" and "approved" files in approval_tests/ sub-directory, keeping the test source directory tidy
auto directoryDisposer = ApprovalTests::Approvals::useApprovalsSubdirectory("approval_tests");
