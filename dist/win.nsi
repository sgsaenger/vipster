Name "Vipster"
!define INSTALLATIONNAME "Vipster"
InstallDir $PROGRAMFILES\${INSTALLATIONNAME}

Page license
Page components
Page directory
Page instfiles
UninstPage uninstConfirm
UninstPage instfiles

LicenseData "license.txt"

Section ""
    SetOutPath $INSTDIR
    File "${buildDirectory}\*"
    WriteUninstaller $INSTDIR\uninstall.exe
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${INSTALLATIONNAME}" "DisplayName" "Vipster"
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${INSTALLATIONNAME}" "UninstallString" '"$INSTDIR\uninstall.exe"'
    WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${INSTALLATIONNAME}" "NoModify" 1
    WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${INSTALLATIONNAME}" "NoRepair" 1
SectionEnd

Section "Start Menu Shortcuts"
    CreateDirectory "$SMPROGRAMS\${INSTALLATIONNAME}"
    CreateShortcut "$SMPROGRAMS\${INSTALLATIONNAME}\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0
    CreateShortcut "$SMPROGRAMS\${INSTALLATIONNAME}\Vipster.lnk" "$INSTDIR\vipster.exe" "" "$INSTDIR\vipster.exe" 0
SectionEnd

Section "Uninstall"
    DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${INSTALLATIONNAME}"
    Delete $INSTDIR\*
    RMDir $INSTDIR
    Delete $SMPROGRAMS\${INSTALLATIONNAME}\*
    RMDir $SMPROGRAMS\${INSTALLATIONNAME}
SectionEnd
