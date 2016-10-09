import os
import sipconfig

# Get the SIP configuration information.
config = sipconfig.Configuration()

# Run SIP to generate the code.
os.system(" ".join([config.sip_bin, "-c", ".", "-b", "vipster.sbf", "vipster.sip"]))

# Create the Makefile.
makefile = sipconfig.SIPModuleMakefile(config, "vipster.sbf")

# Add the library we are wrapping.  The name doesn't include any platform
# specific prefixes or extensions (e.g. the "lib" prefix on UNIX, or the
# ".dll" extension on Windows).
makefile.extra_libs = ["vipster"]

# Generate the Makefile itself.
makefile.generate()
