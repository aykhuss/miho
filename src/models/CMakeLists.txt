add_library(models "")

target_sources(models
  PRIVATE
    Util.cpp
    Model.cpp
    GeometricModel.cpp
    ABCModel.cpp
    ABCNumericModel.cpp
    Scale1DGeometricModel.cpp
    Scale2DGeometricModel.cpp
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/Util.h
    ${CMAKE_CURRENT_LIST_DIR}/Model.h
    ${CMAKE_CURRENT_LIST_DIR}/GeometricModel.h
    ${CMAKE_CURRENT_LIST_DIR}/ABCModel.h
    ${CMAKE_CURRENT_LIST_DIR}/ABCNumericModel.h
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
  cubature
  # adsj
)
