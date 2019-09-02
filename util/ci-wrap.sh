#!/usr/bin/bash

case $1 in
    install)
        echo Installing requirements
        case $2 in
            linux)
                case $3 in
                    web)
                        # setup emscripten
                        git clone --depth 1 --branch incoming https://github.com/urho3d/emscripten-sdk.git ~/emscripten-sdk
                        ~/emscripten-sdk/emsdk activate --build=Release sdk-incoming-64bit binaryen-master-64bit
                        source ~/emscripten-sdk/emsdk_env.sh
                        for compiler in $EMSCRIPTEN/{emcc,em++}; do touch -d "2017-01-01 00:00:00 +0800" $compiler; done
                        ;;
                    desktop)
                        # select python version:
                        pyenv shell 3.7
                        pyenv versions
                        # make sure GCC7 is used:
                        sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 60 --slave /usr/bin/g++ g++ /usr/bin/g++-7
                        sudo update-alternatives --set gcc /usr/bin/gcc-7
                        export PATH="/opt/qt512/bin":$PATH
                        ;;
                esac
                ;;
            osx)
                # install qt
                brew update
                brew install grep
                brew install qt
                export PATH=/usr/local/qt/bin:$PATH
                ;;
            windows)
                # install qt+mingw
                wget "http://download.qt.io/official_releases/online_installers/qt-unified-windows-x86-online.exe" -O qt.exe
                ./qt.exe --verbose --script util/qt-headless.qs
                export MWDIR="/c/Users/travis/Qt/Tools/mingw730_64"
                export QTDIR="/c/Users/travis/Qt/5.12.3/mingw73_64"
                export PATH="$MWDIR/bin:$QTDIR/bin:$PATH"
                # install Python 3.7
                choco install python
                export PATH="/c/Python37:/c/Python37/Scripts:$PATH"
                ;;
        esac
        ;;
    script)
        echo Building and testing
        case $2 in
            linux)
                case $3 in
                    web)
                        mkdir build
                        cd build
                        emcmake cmake -D WEB=YES -D CMAKE_BUILD_TYPE=Release $TRAVIS_BUILD_DIR
                        make -j2
                        ;;
                    desktop)
                        mkdir build
                        cd build
                        cmake -D DESKTOP=YES -D PYTHON=YES -D TESTS=YES -D CMAKE_BUILD_TYPE=Debug -D CMAKE_CXX_FLAGS="-g -O0 -fprofile-arcs -ftest-coverage" $TRAVIS_BUILD_DIR
                        make -j2
                        ./test_lib
                        ;;
                esac
                ;;
            osx)
                mkdir build
                cd build
                cmake -D CMAKE_PREFIX_PATH=/usr/local/opt/qt -D TESTS=YES -D DESKTOP=YES -D CMAKE_BUILD_TYPE=Release $TRAVIS_BUILD_DIR
                make -j2
                ./test_lib
                ;;
            windows)
                mkdir build
                cd build
                mv "C:\Program Files\Git\usr\bin\sh.exe" "C:\Program Files\Git\usr\bin\sh2.exe"
                cmake -D TESTS=YES -D PYTHON=YES -D DESKTOP=YES -D CMAKE_PREFIX_PATH="$QTDIR;$MWDIR" -G "MinGW Makefiles" -D CMAKE_BUILD_TYPE=RELEASE ..
                cmake --build .
                ./test_lib.exe
                ;;
        esac
        ;;
    after_success)
        case $2 in
            linux)
                if [[ $3 = "desktop" ]]
                then
                    echo Collecting coverage
                    bash <(curl -s https://codecov.io/bash) -R $TRAVIS_BUILD_DIR -x gcov-7
                    if [[ $TRAVIS_BRANCH == "testing" ]]
                    then
                        # build AppImage
                        bash $TRAVIS_BUILD_DIR/util/make-appimage.sh
                        # upload continuous build
                        wget https://github.com/d1vanov/ciuploadtool/releases/download/continuous-master/ciuploadtool_linux.zip
                        unzip ciuploadtool_linux.zip
                        chmod 755 ciuploadtool
                        ./ciuploadtool $TRAVIS_BUILD_DIR/Vipster-Linux-x86_64.AppImage
                    fi
                fi
                ;;
            osx)
                if [[ $TRAVIS_BRANCH == "testing" ]]
                then
                    # prepare .dmg file
                    bash $TRAVIS_BUILD_DIR/util/make-osxapp.sh
                    # upload continuous build
                    wget https://github.com/d1vanov/ciuploadtool/releases/download/continuous-master/ciuploadtool_mac.zip
                    unzip ciuploadtool_mac.zip
                    chmod 755 ciuploadtool
                    ./ciuploadtool $TRAVIS_BUILD_DIR/build/Vipster-OSX.dmg
                fi
                ;;
        esac
        ;;
    before_deploy)
        echo Preparing deployment
        case $2 in
            linux)
                case $3 in
                    web)
                        # prepare website
                        mv vipster.{js,wasm} $TRAVIS_BUILD_DIR/gh-pages/emscripten
                        ;;
                    desktop)
                        # build AppImage
                        bash $TRAVIS_BUILD_DIR/util/make-appimage.sh
                        # enable deployment
                        export DEPLOY_FILE=$TRAVIS_BUILD_DIR/Vipster-Linux-x86_64.AppImage
                        ;;
                esac
                ;;
            osx)
                # prepare .dmg file
                bash $TRAVIS_BUILD_DIR/util/make-osxapp.sh
                # enable deployment
                export DEPLOY_FILE=$TRAVIS_BUILD_DIR/build/Vipster-OSX.dmg
                ;;
        esac
        ;;
esac
