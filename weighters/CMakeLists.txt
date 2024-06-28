file(GLOB HEADER_FILES *.h)

install(FILES ${HEADER_FILES} DESTINATION include/weighters)

install(FILES scale_factors.root DESTINATION etc)
