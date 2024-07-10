set (
    SRCS
    fips202.cpp
    poly.cpp
    polyvec.cpp
    polymat.cpp
    sample.cpp
)

set(
    HEADERS
    fips202.h
    poly.h
    polymat.h
    sample.h
)

add_library (
    poly
    STATIC
    ${SRCS}
    ${HEADERS}
)

