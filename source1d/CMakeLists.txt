SET(APP_NAME source1d)

# TODO - Add an entry here for every source file in your project. Must not be empty!
SET(SourceFileNames
    source1d.cpp
    main.cpp

startUp1d.cpp
)

# TODO - Add an entry here for every header file in your project. Must not be empty!
SET(HeaderFileNames
source1d.h

startUp1d.h
all.h

)

ADD_EXECUTABLE(${APP_NAME} ${SourceFileNames} ${HeaderFileNames})


