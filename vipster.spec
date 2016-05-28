# -*- mode: python -*-

block_cipher = None


a = Analysis(['scripts\\vipster'],
             pathex=['C:\\Users\\Hein\\Desktop\\vipster'],
             binaries=None,
             datas=[('vipster/default.json','vipster'),
                    ('vipster/gui/opengl/*.vert','vipster/gui/opengl'),
                    ('vipster/gui/opengl/*.frag','vipster/gui/opengl'),
                    ('vipster/gui/opengl/sphere_model','vipster/gui/opengl'),
                    ('vipster/gui/opengl/bond_model','vipster/gui/opengl'),
                    ('README.md','.'),('LICENSE','.')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
#          exclude_binaries=True,
          name='vipster',
          debug=False,
          strip=False,
          upx=True,
          console=True , icon='vipster_icon.ico')
#coll = COLLECT(exe,
#               a.binaries,
#               a.zipfiles,
#               a.datas,
#               strip=False,
#               upx=True,
#               name='vipster')
