add_library(models "")

# find_package(Boost)
# if(Boost_FOUND)
#   include_directories(${Boost_INCLUDE_DIRS})
# endif()

target_sources(models
  PRIVATE
    Util.cpp
    Model.cpp
    GeometricModel.cpp
    ABCModel.cpp
    ABCNumericModel.cpp
    Scale1DModel.cpp
    Scale2DModel.cpp
    Scale1DGeometricModel.cpp
    Scale2DGeometricModel.cpp
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/Util.h
    ${CMAKE_CURRENT_LIST_DIR}/Model.h
    ${CMAKE_CURRENT_LIST_DIR}/GeometricModel.h
    ${CMAKE_CURRENT_LIST_DIR}/ABCModel.h
    ${CMAKE_CURRENT_LIST_DIR}/ABCNumericModel.h
    ${CMAKE_CURRENT_LIST_DIR}/Scale1DModel.h
    ${CMAKE_CURRENT_LIST_DIR}/Scale2DModel.h
    ${CMAKE_CURRENT_LIST_DIR}/Scale1DGeometricModel.h
    ${CMAKE_CURRENT_LIST_DIR}/Scale2DGeometricModel.h
  )
target_include_directories(models
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )
target_link_libraries(models
  PUBLIC
  fmt::fmt
  PRIVATE
  cubature
  # Boost::headers
  # adsj
)
