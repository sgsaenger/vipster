#!/usr/bin/bash

case $1 in
    install)
        echo Installing requirements
        case $2 in
            linux)
                case $3 in
                    web)
                        cd $HOME
                        # setup emscripten
                        git clone --depth 1 https://github.com/emscripten-core/emsdk.git
                        $HOME/emsdk/emsdk install latest
                        $HOME/emsdk/emsdk activate latest
                        source $HOME/emsdk/emsdk_env.sh
                        ;;
                    desktop)
                        cd $HOME
                        # install qt
                        pip3 install aqtinstall
                        aqt install -O $HOME/Qt 5.14.0 linux desktop gcc_64
                        export QTDIR="$HOME/Qt/5.14.0/gcc_64"
                        export PATH="$QTDIR/bin:$PATH"
                        ;;
                esac
                ;;
            osx)
                cd $HOME
                # install qt
                pip3 install aqtinstall
                aqt install -O $HOME/Qt 5.14.0 mac desktop clang_64
                export QTDIR="$HOME/Qt/5.14.0/clang_64"
                export PATH="$QTDIR/bin:$PATH"
                ;;
            windows)
                cd $HOME
                # install nuwen-mingw
                wget https://nuwen.net/files/mingw/mingw-17.0-without-git.exe
                7z x mingw-17.0-without-git.exe
                MinGW/set_distro_paths.bat
                export MWDIR="$HOME/MinGW"
                # install qt
                pip3 install aqtinstall
                aqt install -O $HOME/Qt 5.14.0 windows desktop win64_mingw73
                export QTDIR="$HOME/Qt/5.14.0/mingw73_64"
                export PATH="$MWDIR/bin:$QTDIR/bin:$PATH"
                ;;
        esac
        ;;
    script)
        echo Building and testing
        case $2 in
            linux)
                case $3 in
                    web)
                        mkdir -p $HOME/build/release
                        cd $HOME/build/release
                        emcmake cmake -DWEB=ON -DCMAKE_BUILD_TYPE=Release $SOURCE_DIR
                        make -j2
                        ;;
                    desktop)
                        mkdir -p $HOME/build/debug
                        cd $HOME/build/debug
                        cmake -DDESKTOP=ON -DPYSHELL=ON -DPYBIND=ON -DTESTS=ON -DCMAKE_BUILD_TYPE=Debug -DPYTHON_EXECUTABLE=$PY_DIR/bin/python3 -DCMAKE_CXX_FLAGS="-g -O0 -fprofile-arcs -ftest-coverage" $SOURCE_DIR
                        make -j2
                        ./test_lib
                        ;;
                esac
                ;;
            osx)
                mkdir -p $HOME/build/release
                cd $HOME/build/release
                cmake -DCMAKE_PREFIX_PATH=$QTDIR -DTESTS=ON -DDESKTOP=ON -DCMAKE_BUILD_TYPE=Release $SOURCE_DIR
                make -j2
                ./test_lib
                ;;
            windows)
                mkdir -p $HOME/build/release
                cd $HOME/build/release
                mv "C:\Program Files\Git\usr\bin\sh.exe" "C:\Program Files\Git\usr\bin\sh2.exe"
                cmake -DTESTS=ON -DDESKTOP=ON -DPYSHELL=ON -DPYBIND=ON -DCMAKE_PREFIX_PATH="$QTDIR;$MWDIR" -G"MinGW Makefiles" -DCMAKE_BUILD_TYPE=RELEASE ..
                cmake --build . --
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
                    bash <(curl -s https://codecov.io/bash) -R $SOURCE_DIR -x gcov-8
                    cd $SOURCE_DIR
                    BRANCH=$(git symbolic-ref --short -q HEAD)
                    if [[ $BRANCH == "CI" ]]
                    then
                        # build release-version
                        mkdir -p $HOME/build/release
                        cd $HOME/build/release
                        cmake -DDESKTOP=ON -DPYSHELL=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DPYTHON_EXECUTABLE=$PY_DIR/bin/python3 $SOURCE_DIR
                        make -j2
                        # build AppImage
                        bash $SOURCE_DIR/util/make-appimage.sh
                        # upload continuous build
                        wget https://github.com/d1vanov/ciuploadtool/releases/download/continuous-master/ciuploadtool_linux.zip
                        unzip ciuploadtool_linux.zip
                        chmod 755 ciuploadtool
                        ./ciuploadtool $SOURCE_DIR/release/Vipster-Linux-x86_64.AppImage -suffix test
                    fi
                fi
                ;;
            osx)
                cd $SOURCE_DIR
                BRANCH=$(git symbolic-ref --short -q HEAD)
                if [[ $BRANCH == "CI" ]]
                then
                    # prepare .dmg file
                    cd $SOURCE_DIR/build
                    bash $SOURCE_DIR/util/make-osxapp.sh
                    # upload continuous build
                    wget https://github.com/d1vanov/ciuploadtool/releases/download/continuous-master/ciuploadtool_mac.zip
                    unzip ciuploadtool_mac.zip
                    chmod 755 ciuploadtool
                    ./ciuploadtool $SOURCE_DIR/build/Vipster-OSX.dmg -suffix test
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
                        mv vipster.{js,wasm} $SOURCE_DIR/gh-pages/emscripten
                        ;;
                    desktop)
                        # build release-version
                        mkdir -p $SOURCE_DIR/release
                        cd $SOURCE_DIR/release
                        cmake -DDESKTOP=ON -DPYSHELL=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr $SOURCE_DIR
                        make -j2
                        # build AppImage
                        bash $SOURCE_DIR/util/make-appimage.sh
                        # enable deployment
                        export DEPLOY_FILE=$SOURCE_DIR/release/Vipster-Linux-x86_64.AppImage
                        ;;
                esac
                ;;
            osx)
                # prepare .dmg file
                bash $SOURCE_DIR/util/make-osxapp.sh
                # enable deployment
                export DEPLOY_FILE=$SOURCE_DIR/build/Vipster-OSX.dmg
                ;;
        esac
        ;;
esac
