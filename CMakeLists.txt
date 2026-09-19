cmake_minimum_required(VERSION 3.6.0)
project(zp_physics)


include_directories(
	${PROJECT_SOURCE_DIR}/include
	${PROJECT_SOURCE_DIR}/src
)


file(GLOB SRCS
	"${PROJECT_SOURCE_DIR}/include/*"
	"${PROJECT_SOURCE_DIR}/src/*"
	"${PROJECT_SOURCE_DIR}/src/zp_physics/*"
	"${PROJECT_SOURCE_DIR}/src/zp_physics/core2d/*"
)


add_library(zp_physics SHARED ${SRCS})
add_library(zp_physics STATIC ${SRCS})

