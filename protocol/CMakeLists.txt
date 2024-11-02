
set (
    funcSRC
    main_protocol.cpp
    pack.cpp
    protocol1.cpp
    protocol2.cpp
    cpaenc.cpp
    tools.cpp
    construct.cpp
)
set (
    mainSRC
    app.cpp
    ${funcSRC}
)
set (
    multiSRC
    multi_protocol.cpp
    protocol_batch.cpp
    pack.cpp
    protocol1.cpp
    cpaenc.cpp
    tools.cpp
    construct.cpp
)
set (
    multitestSRC
    speed_test_batch.cpp
    protocol_batch.cpp
    pack.cpp
    protocol1.cpp
    cpaenc.cpp
    tools.cpp
    construct.cpp
)
set(
    speedtestSRC
    speed_test.cpp
    ${funcSRC}
)
add_executable(main_protocol ${mainSRC})
target_link_libraries(main_protocol PRIVATE poly)
add_executable(speed_test ${speedtestSRC})
target_link_libraries(speed_test PRIVATE poly)
add_executable(multi_protocol ${multiSRC})
target_link_libraries(multi_protocol PRIVATE poly)

add_executable(speed_test_batch ${multitestSRC})
target_link_libraries(speed_test_batch PRIVATE poly)