TEMPLATE = subdirs

SUBDIRS += \
    libvipster \
    vipster \
    tests

unix: SUBDIRS += python
