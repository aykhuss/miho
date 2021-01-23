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
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/Util.h
    ${CMAKE_CURRENT_LIST_DIR}/Model.h
    ${CMAKE_CURRENT_LIST_DIR}/ModelPrototype.h
    ${CMAKE_CURRENT_LIST_DIR}/GeometricModel.h
    ${CMAKE_CURRENT_LIST_DIR}/ABCModel.h
    ${CMAKE_CURRENT_LIST_DIR}/ABCNumericModel.h
    ${CMAKE_CURRENT_LIST_DIR}/ScaleModel.h
    ${CMAKE_CURRENT_LIST_DIR}/ScaleModel.hpp
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
