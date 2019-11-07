# add_subdirectory(Havec)
file(COPY Havec/Havec DESTINATION ${CMAKE_BINARY_DIR}/GraphConstruction/Havec)

#add_subdirectory(bcalm)
file(COPY bcalm/build/bcalm DESTINATION ${CMAKE_BINARY_DIR}/GraphConstruction/bcalm)