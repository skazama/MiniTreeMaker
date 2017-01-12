project(MyAnalysis)

add_library(mytool AnalysisTool.cxx) 
add_library(myana MyAnalysis.cxx) 
add_library(json jsoncpp.cpp) 
add_library(paxjson GetJSON.cxx) 
add_library(dict Dict.cc) 
add_executable(ana main.cxx) 

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${MY_ANALYSIS_INCLUDE_DIRS})
link_directories(${CMAKE_SOURCE_DIR}/build/lib/)

target_link_libraries(ana LINK_PUBLIC ${ROOT_LIBRARIES} ${MPI_LIBRARIES} /home/atp/kazama/xenon1t/MiniTreeMaker/v0.7/include/pax.so)
target_link_libraries(ana LINK_PUBLIC ${ROOT_LIBRARIES} ${MPI_LIBRARIES} myana)
target_link_libraries(ana LINK_PUBLIC ${ROOT_LIBRARIES} ${MPI_LIBRARIES} mytool)
target_link_libraries(ana LINK_PUBLIC ${ROOT_LIBRARIES} ${MPI_LIBRARIES} json)
target_link_libraries(ana LINK_PUBLIC ${ROOT_LIBRARIES} ${MPI_LIBRARIES} paxjson)
target_link_libraries(ana LINK_PUBLIC ${ROOT_LIBRARIES} ${MPI_LIBRARIES} dict)
