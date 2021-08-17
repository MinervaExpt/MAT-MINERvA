file(GLOB SOURCE_FILES *.cxx)
file(GLOB HEADER_FILES *.h)

install(FILES ${SOURCE_FILES} ${HEADER_FILES} DESTINATION include/utilities)
