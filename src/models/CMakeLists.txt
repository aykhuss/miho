add_library(models "")

target_sources(models
  PRIVATE
    Model.cpp
    GeometricModel.cpp
    ABCModel.cpp
    Scale1DGeometricModel.cpp
    Scale2DGeometricModel.cpp
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/Model.h
    ${CMAKE_CURRENT_LIST_DIR}/GeometricModel.h
    ${CMAKE_CURRENT_LIST_DIR}/ABCModel.h
    ${CMAKE_CURRENT_LIST_DIR}/Scale1DGeometricModel.h
    ${CMAKE_CURRENT_LIST_DIR}/Scale2DGeometricModel.h
  )
target_include_directories(models
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )
target_link_libraries(models
	PRIVATE
	fmt::fmt
)
