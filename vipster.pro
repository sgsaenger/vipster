TEMPLATE = subdirs

SUBDIRS += \
    libvipster \
    vipster
#unix:  SUBDIRS += python
!wasm: SUBDIRS += tests
