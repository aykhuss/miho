set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

add_executable(miho main.cpp)

add_subdirectory(models)

target_link_libraries(miho
  PUBLIC
  fmt::fmt
	PRIVATE
	models
  CLI11::CLI11
)

