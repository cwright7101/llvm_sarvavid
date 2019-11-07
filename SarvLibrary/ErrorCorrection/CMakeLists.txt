#Rcorrector
file(COPY rcorrector DESTINATION ${CMAKE_BINARY_DIR}/ErrorCorrection/)

#bloocoo
file(COPY bloocoo/build/bin DESTINATION ${CMAKE_BINARY_DIR}/ErrorCorrection/bloocoo/)

#coral
file(COPY coral/coral DESTINATION ${CMAKE_BINARY_DIR}/ErrorCorrection/coral/)

#LoRDEC
file(COPY LoRDEC/build/tools DESTINATION ${CMAKE_BINARY_DIR}/ErrorCorrection/LoRDEC/)