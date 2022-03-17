#include <emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/html5.h>

#include "guiwrapper.h"
#include "vipster/molecule.h"
#include "vipster/fileio.h"
#include "vipster/plugins/json.h"

namespace em = emscripten;
using namespace Vipster;

static PluginList plugins = defaultPlugins();

EM_BOOL mouse_event(int eventType, const EmscriptenMouseEvent* mouseEvent, void *gui_ptr)
{
    if(mouseEvent->button == 2){
        return 0; // don't act on right click
    }
    auto &gui = *reinterpret_cast<GuiWrapper*>(gui_ptr);
    enum class MouseMode { Camera, Select, Modify};
    enum class OpMode { None, Rotation, Translation };
    static OpMode currentOp = OpMode::None;
    static long localX, localY;
    switch (eventType) {
    case EMSCRIPTEN_EVENT_MOUSEDOWN:
        if(currentOp == OpMode::None){
            int button = mouseEvent->button | mouseEvent->shiftKey | mouseEvent->ctrlKey << 1;
            switch(button){
            case 0:
                currentOp = OpMode::Rotation;
                break;
            case 1:
                currentOp = OpMode::Translation;
                break;
            case 2:
                gui.alignViewMat(GuiWrapper::alignDir::z);
            }
        }
        break;
    case EMSCRIPTEN_EVENT_MOUSEUP:
        if(currentOp!=OpMode::None && !mouseEvent->buttons){
            currentOp = OpMode::None;
        }
        break;
    case EMSCRIPTEN_EVENT_MOUSELEAVE:
        currentOp = OpMode::None;
        break;
    case EMSCRIPTEN_EVENT_MOUSEMOVE:
        switch(currentOp){
        case OpMode::Rotation:
            gui.rotateViewMat(mouseEvent->clientX-localX,
                              mouseEvent->clientY-localY, 0);
            break;
        case OpMode::Translation:
            gui.translateViewMat(mouseEvent->clientX-localX,
                                 -(mouseEvent->clientY-localY), 0);
            break;
        default:
            break;
        }
        break;
    }
    localX = mouseEvent->clientX;
    localY = mouseEvent->clientY;
    return 1;
}

EM_BOOL wheel_event(int, const EmscriptenWheelEvent* wheelEvent, void *gui_ptr)
{
    auto &gui = *reinterpret_cast<GuiWrapper*>(gui_ptr);
    gui.zoomViewMat(wheelEvent->deltaY<0?1.1:0.9);
    return 1;
}

EM_BOOL draw_loop(double time, void *view_ptr);

class VipsterView : public GuiWrapper
{
public:
    explicit VipsterView(const std::string& canvasID)
    : GuiWrapper{Vipster::settings},
      canvas{canvasID}
    {
        // initialize GL context
        EmscriptenWebGLContextAttributes attrs;
        emscripten_webgl_init_context_attributes(&attrs);
        attrs.enableExtensionsByDefault = 1;
        attrs.preserveDrawingBuffer = true;
        attrs.majorVersion = 2;
        attrs.minorVersion = 0;
        context = emscripten_webgl_create_context(canvas.c_str(), &attrs);
        if (context<=0)
        {
            printf("Could not create WebGL2 context on canvas %s: %d\n", canvas.c_str(), context);
            throw Error("Could not create WebGL2 context on canvas " + canvas + ": " + std::to_string(context));
        }else{
            printf("Creating context on canvas %s\n", canvas.c_str());
        }
        emscripten_webgl_make_context_current(context);
        initGL();

        // register mouse handling
        emscripten_set_mousedown_callback(canvas.c_str(), this, 0, mouse_event);
        emscripten_set_mouseup_callback(canvas.c_str(), this, 0, mouse_event);
        emscripten_set_mousemove_callback(canvas.c_str(), this, 0, mouse_event);
        emscripten_set_mouseleave_callback(canvas.c_str(), this, 0, mouse_event);
        emscripten_set_wheel_callback(canvas.c_str(), this, 0, wheel_event);

        // register simple draw-loop (can exist multiple times)
        emscripten_request_animation_frame_loop(draw_loop, this);
    }
    ~VipsterView()
    {
        emscripten_webgl_destroy_context(context);
    }
    void draw()
    {
        if(!curStep) return;
        emscripten_webgl_make_context_current(context);
        try{
            // handle resize
            int canvas_width{}, canvas_height{};
            double css_width{}, css_height{};
            emscripten_get_canvas_element_size(canvas.c_str(), &canvas_width, &canvas_height);
            emscripten_get_element_css_size(canvas.c_str(), &css_width, &css_height);
            if(canvas_width != css_width || canvas_height != css_height){
                canvas_width = css_width;
                canvas_height = css_height;
                emscripten_set_canvas_element_size(canvas.c_str(), canvas_width, canvas_height);
            }
            // actual draw
            emscripten_webgl_make_context_current(context);
            resizeViewMat(canvas_width, canvas_height);
            GuiWrapper::draw(reinterpret_cast<void*>(context));
        }
        catch(Vipster::Error e){
            printf("%s\n", e.what());
        }
        catch(const std::exception& e){
            printf("%s\n", e.what());
        }
        catch(...){
            printf("Unknown error\n");
        }
    }
private:
    std::string canvas;
    EMSCRIPTEN_WEBGL_CONTEXT_HANDLE context;
};

EM_BOOL draw_loop(double time, void *view_ptr)
{
    auto &view = *reinterpret_cast<VipsterView*>(view_ptr);
    view.draw();
    return EM_TRUE;
}

// wrap to IO stuff
size_t nplug(){ return plugins.size(); }
std::string plugName(int i){ return plugins[i]->name; }
bool plugRead(int i){ return static_cast<bool>(plugins[i]->parser); }
bool plugWrite(int i){ return static_cast<bool>(plugins[i]->writer); }
int guessFmtEm(std::string s){
    auto pos = s.find_last_of('.');
    if(pos != s.npos){
        s = s.substr(pos+1);
    }else{
        pos = s.find_last_of("/\\");
        if(pos != s.npos){
            s = s.substr(pos+1);
        }
    }
    auto plug = std::find_if(plugins.begin(), plugins.end(), [&](const Plugin *p){
        return p->extension == s;
    });
    if(plug != plugins.end()){
        return plug-plugins.begin();
    }else{
        return -1;
    }
}
std::string addExtension(std::string s, int f){
    auto ext = plugins[f]->extension;
    auto pos = s.find_last_of('.');
    if((pos == s.npos) && (s.substr(pos+1) != ext)){
        return s.substr(0, pos) + '.' + ext;
    }
    return s;
}
std::string emWriteFile(const Molecule &m, int s, int f){
    try{
        const auto& plug = plugins[f];
        Parameter param{nullptr};
        if(plug->makeParam){
            param = plug->makeParam();
        }
        Preset preset{nullptr};
        if(plug->makePreset){
            preset = plug->makePreset();
        }
        writeFile("/tmp/output.file", plug, m,
                  (size_t)s, param, preset);
        return "";
    } catch(std::exception &e) {
        return e.what();
    }
}


EMSCRIPTEN_BINDINGS(vipster)
{
    // expose IO stuff via free functions
    em::function("nplug", &nplug);
    em::function("plugName", &plugName);
    em::function("plugRead", &plugRead);
    em::function("plugWrite", &plugWrite);
    em::function("guessFmt", &guessFmtEm);
    em::function("addExtension", &addExtension);
    em::function("writeFile", &emWriteFile);
    // GUI Wrapper
    em::class_<VipsterView>("VipsterView")
        .constructor<const std::string&>()
        .function("setStep", std::function([](VipsterView &v, Step& s){ v.setMainStep(&s); }))
        .function("draw", &VipsterView::draw)
        .property("mult",
                  std::function([](const VipsterView &v){ return v.mult; }),
                  std::function([](VipsterView &v, GUI::PBCVec m){ v.mult = m; }))
        .function("zoom", std::function([](VipsterView &v, float val){ v.zoomViewMat(val); }))
        .function("rotate", std::function([](VipsterView &v, int x, int y){ v.rotateViewMat(x, y, 0); }))
        .function("translate", std::function([](VipsterView &v, int x, int y){ v.translateViewMat(x, y, 0); }))
        ;
    // libvipster API
    em::class_<Molecule>("Molecule")
        .constructor(std::function([](){ return Molecule(); }))
        // create from JSON-string
        .constructor(std::function([](std::string s) -> Molecule{
            auto ss = std::istringstream(s);
            return std::get<0>(Plugins::JSON.parser("", ss));
        }))
        // create from File
        .constructor(std::function([](const std::string& fn, int fmt){
            return std::get<0>(readFile(fn, plugins[fmt]));
        }))
        .property("name", &Molecule::name)
        .function("getStep",
                  std::function([](Molecule *m, size_t idx){ return &m->getStep(idx); }),
                  em::allow_raw_pointers()) // need to return a pointer so embind won't create a copy of Step
        .property("nstep", std::function([](const Molecule& m){return m.getNstep();}))
        ;
    em::class_<Step>("Step")
        .property("nat", std::function([](const Step& s){ return s.getNat(); }))
        .function("getAtom", std::function([](Step &s, size_t i){ return s.at(i); }))
        .function("getAtomIt", std::function([](Step &s){ return s.begin(); }))
        .property("fmt",
                  std::function([](const Step& s){ return static_cast<int>(s.getFmt());}),
                  std::function([](Step &s, int f){ s.setFmt(static_cast<AtomFmt>(f), true);}))
        .property("hasCell", &Step::hasCell, std::function([](Step& s, bool b){s.enableCell(b);}))
        .property("cellDim",
                  std::function([](const Step& s){ return s.getCellDim(AtomFmt::Angstrom); }),
                  std::function([](Step &s, double d){ s.setCellDim(d, AtomFmt::Angstrom, false); }))
        .property("cellVec",
                  std::function([](const Step& s){ return s.getCellVec(); }),
                  std::function([](Step &s, Mat m){ s.setCellVec(m, false); }))
        .function("generateBonds", std::function([](Step &s){ s.generateBonds(); }))
        .function("hasBonds", std::function([](Step &s){ return !s.getBonds().empty(); }))
        ;
    em::class_<Step::atom>("Atom")
        .property("name",
                  std::function([](const Step::atom& at)->std::string{ return at.name; }),
                  std::function([](Step::atom& at, const std::string &name){ at.name = name; }))
        .property("coord",
                  std::function([](const Step::atom& at)->Vec{ return at.coord; }),
                  std::function([](Step::atom& at, const Vec &coord){ at.coord = coord; }))
        ;
    em::class_<Step::iterator>("AtomIterator")
        .function("increment", em::select_overload<Step::iterator&()>(&Step::iterator::operator++))
        .property("name",
                  std::function([](const Step::iterator& it)->std::string{ return it->name; }),
                  std::function([](Step::iterator& it, const std::string &name){ it->name = name; }))
        .property("coord",
                  std::function([](const Step::iterator& it)->Vec{ return it->coord; }),
                  std::function([](Step::iterator& it, const Vec &coord){ it->coord = coord; }))
        ;
    em::value_array<GUI::PBCVec>("PBCVec")
            .element(em::index<0>())
            .element(em::index<1>())
            .element(em::index<2>());
    em::value_array<Vec>("Vec")
            .element(em::index<0>())
            .element(em::index<1>())
            .element(em::index<2>());
    em::value_array<Mat>("Mat")
            .element(em::index<0>())
            .element(em::index<1>())
            .element(em::index<2>());
}
