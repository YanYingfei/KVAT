
set (
    mainSRC
    main_protocol.cpp
    pack.cpp
    protocol1.cpp
    protocol2.cpp
    cpaenc.cpp
    tools.cpp
    construct.cpp
)


add_executable(main_protocol ${mainSRC})
target_link_libraries(main_protocol PRIVATE poly)