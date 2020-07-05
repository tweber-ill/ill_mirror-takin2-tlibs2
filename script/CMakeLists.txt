#
# @author Tobias Weber
# @date 11-apr-2020
# @license: see 'LICENSE.GPL' file
#

project(mcalc)
cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE ON)


find_package(FLEX 2 REQUIRED)
find_package(BISON 3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)


add_definitions(-std=c++17)
add_definitions(${Boost_CXX_FLAGS})

include_directories(
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_SOURCE_DIR}/src"
	"${PROJECT_SOURCE_DIR}/.."
	"${FLEX_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/.."
	"${CMAKE_CURRENT_BINARY_DIR}"
)



BISON_TARGET(parser_impl
	src/parser.y ${CMAKE_CURRENT_BINARY_DIR}/parser_impl.cpp
	DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/parser_defs.h
	VERBOSE REPORT_FILE ${CMAKE_CURRENT_BINARY_DIR}/parser_report.txt
)


FLEX_TARGET(lexer_impl
	src/lexer.l ${CMAKE_CURRENT_BINARY_DIR}/lexer_impl.cpp
	DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/lexer_impl.h
)

ADD_FLEX_BISON_DEPENDENCY(lexer_impl parser_impl)



message("Flex outputs: ${FLEX_lexer_impl_OUTPUTS}")
message("Bison output source: ${BISON_parser_impl_OUTPUT_SOURCE}")
message("Bison output header: ${BISON_parser_impl_OUTPUT_HEADER}")



add_executable(mcalc
	src/parser.cpp src/parser.h src/ast.h src/sym.h
	src/llasm.cpp src/llasm.h
	src/printast.cpp src/printast.h
	${FLEX_lexer_impl_OUTPUTS}
	${BISON_parser_impl_OUTPUT_SOURCE} ${BISON_parser_impl_OUTPUT_HEADER}
	../libs/log.cpp ../libs/log.h
)

#add_dependencies(parser parser_impl lexer_impl)
target_link_libraries(mcalc ${Boost_LIBRARIES})



# build runtime
execute_process(COMMAND clang -O2 -S -emit-llvm -o runtime_opt.asm "${PROJECT_SOURCE_DIR}/src/runtime.cpp")
execute_process(COMMAND clang -S -emit-llvm -o runtime.asm "${PROJECT_SOURCE_DIR}/src/runtime.cpp")
