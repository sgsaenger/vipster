from imp import reload


def test_userconf():
    from vipster import settings
    settings.saveConfig()
    reload(settings)
