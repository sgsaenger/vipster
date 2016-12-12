#include "ioplugins/pwinput.h"
#define BOOST_SPIRIT_DEBUG
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/fusion/include/std_pair.hpp>

namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;
using namespace std;

std::ostream& operator<< (ostream& stream, const map<string, string>m)
{
    for(auto& p:m){
        stream << p.first << " has value " << p.second << endl;
    }
    return stream;
}

std::ostream& operator<< (ostream& stream, const Vipster::IO::PWParam& p)
{
    stream << "control" << endl << p.control << endl;
    stream << "system" << endl << p.system << endl;
    stream << "electrons" << endl << p.electrons << endl;
    stream << "ions" << endl << p.ions << endl;
    stream << "cell" << endl << p.cell << endl;
    return stream;
}

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::Atom,
        (std::string, name)
        (float, coord[0])
        (float, coord[1])
        (float, coord[2])
)

BOOST_FUSION_ADAPT_ADT(
        Vipster::Step,
        (std::vector<Vipster::Atom>, std::vector<Vipster::Atom>, obj.getAtoms(), obj.newAtoms(val))
)

BOOST_FUSION_ADAPT_ADT(
        Vipster::Molecule,
        (obj.getSteps(), obj.newSteps(val))
        (obj.getKPointFmt(), obj.setKPointFmt(val))
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::PseEntry,
        m,
        PWPP
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::IO::PWParam,
        control,
        system,
        electrons,
        ions,
        cell
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::IO::PWData,
        data,
        mol
)

template<typename Iterator>
struct pwi_parse_grammar
        : qi::grammar<Iterator, Vipster::IO::PWData(), qi::blank_type>
{
    pwi_parse_grammar(): pwi_parse_grammar::base_type(file)
    {
        key = *(qi::graph - qi::char_("=/,'"));
        key.name("Key (or unquoted value)");
        fstr = qi::char_('\'') >> qi::lexeme[*(qi::char_ - '\'')] > qi::char_('\'');
        fstr.name("Quoted value");
        value = fstr | key;
        value.name("Value");
        nl_entry = key >> qi::lit('=') >> value;
        nl_entry.name("Namelist KV-pair");
        namelist = qi::omit[qi::no_case[qi::string(qi::_r1)]] >> *qi::eol
                   >> nl_entry % (*(qi::eol >> -(qi::lit('!') >> *(qi::char_ - qi::eol) > qi::eol)) | ',')
                   >> *qi::eol >> '/';
        namelist.name("Namelist");
        param = namelist(string{"&control"}) > *qi::eol >
                namelist(string{"&system"}) > *qi::eol >
                namelist(string{"&electrons"}) > *qi::eol >
                -(namelist(string{"&ions"}) > *qi::eol) >
                -(namelist(string{"&cell"}));
        param.name("PW parameter set");
        //TODO: hier weitermachen
        //IDEE: semantic action nach param -> nat,ntyp,ibrav an mol uebergeben?!
        //IDEE2: ^-operator um reihenfolge beliebig zu machen?
        name = +(qi::char_ - qi::space);
        pseentry = qi::double_ > name;
        pseentry.name("PSEEntry");
        psepair = name > pseentry;
        psepair.name("PSE KV-pair");
        species = qi::no_case[qi::lit("atomic_species")] > qi::eol
                > (psepair % qi::eol);
        species.name("Atomic species");
        atom = name > qi::double_ > qi::double_ > qi::double_;
        atom.name("Atom");
        format = qi::no_case[qi::lit("crystal")][phx::ref(fmt) = Vipster::AtomFmt::Crystal]
               | qi::no_case[qi::lit("alat")][phx::ref(fmt) = Vipster::AtomFmt::Alat]
               | qi::no_case[qi::lit("bohr")][phx::ref(fmt) = Vipster::AtomFmt::Bohr]
               | qi::no_case[qi::lit("angstrom")][phx::ref(fmt) = Vipster::AtomFmt::Angstrom];
        format.name("Atom format");
        atoms = qi::no_case[qi::lit("atomic_positions")] > -format > qi::eol
                > (atom % qi::eol);
        atoms.name("Atomic positions");
        positions = qi::as<Vipster::Step>()[atoms];
        kpoints %= qi::no_case[qi::lit("k_points")] >
                    (qi::no_case[qi::lit("gamma")][qi::_val = Vipster::KPointFmt::Gamma]
                    |qi::no_case[qi::lit("automatic")][qi::_val = Vipster::KPointFmt::MPG]);
        mol.name("Molecule");
        mol = qi::omit[species] > *qi::eol >
              positions > *qi::eol >
              kpoints > *qi::eol;// >
//              cell > *qi::eol >
//              occupations > *qi::eol >
//              constraints > *qi::eol >
//              forces > *qi::eol;
        file = param > *qi::eol > mol;
        file.name("PWScf input file");
        qi::on_error<qi::fail>(
                file,
                std::cout
                << phx::val("Error! Expecting ")
                << qi::_4
                << phx::val(" here: \"")
                << phx::construct<std::string>(qi::_3, qi::_2)
                << phx::val("\"")
                << std::endl
        );
//        BOOST_SPIRIT_DEBUG_NODE(file);
//        BOOST_SPIRIT_DEBUG_NODE(param);
//        BOOST_SPIRIT_DEBUG_NODE(namelist);
//        BOOST_SPIRIT_DEBUG_NODE(nl_entry);
//        BOOST_SPIRIT_DEBUG_NODE(value);
//        BOOST_SPIRIT_DEBUG_NODE(key);
//        BOOST_SPIRIT_DEBUG_NODE(fstr);
//        BOOST_SPIRIT_DEBUG_NODE(mol);
        BOOST_SPIRIT_DEBUG_NODE(name);
//        BOOST_SPIRIT_DEBUG_NODE(species);
//        BOOST_SPIRIT_DEBUG_NODE(psepair);
//        BOOST_SPIRIT_DEBUG_NODE(pseentry);
        BOOST_SPIRIT_DEBUG_NODE(atoms);
        BOOST_SPIRIT_DEBUG_NODE(atom);
        BOOST_SPIRIT_DEBUG_NODE(format);
    }
    Vipster::AtomFmt fmt = Vipster::AtomFmt::Alat;
    qi::rule<Iterator, Vipster::IO::PWData(), qi::blank_type> file;
//Parameter-set, composed of F90-namelists
    qi::rule<Iterator, Vipster::IO::PWParam(), qi::blank_type> param;
    qi::rule<Iterator, map<string, string>(string), qi::blank_type> namelist;
    qi::rule<Iterator, pair<string, string>(), qi::blank_type> nl_entry;
    qi::rule<Iterator, string()> value;
    qi::rule<Iterator, string()> key;
    qi::rule<Iterator, string()> fstr;
//Molecule
    qi::rule<Iterator, Vipster::Molecule(), qi::blank_type> mol;
    qi::rule<Iterator, vector<Vipster::Step>(), qi::blank_type> positions;
    qi::rule<Iterator, Vipster::Step(), qi::blank_type> step;
    qi::rule<Iterator, string()> name;
    qi::rule<Iterator, vector<pair<string,Vipster::PseEntry>>(), qi::blank_type> species;
    qi::rule<Iterator, pair<string,Vipster::PseEntry>(), qi::blank_type> psepair;
    qi::rule<Iterator, Vipster::PseEntry(), qi::blank_type> pseentry;
    qi::rule<Iterator, vector<Vipster::Atom>(), qi::blank_type> atoms;
    qi::rule<Iterator, qi::blank_type> format;
    qi::rule<Iterator, Vipster::Atom(), qi::blank_type> atom;
    qi::rule<Iterator, Vipster::KPointFmt(), qi::blank_type> kpoints;
//    qi::rule<Iterator, qi::unused_type(), qi::blank_type> cell;
//    qi::rule<Iterator, qi::unused_type(), qi::blank_type> occupations;
//    qi::rule<Iterator, qi::unused_type(), qi::blank_type> constraints;
//    qi::rule<Iterator, qi::unused_type(), qi::blank_type> forces;
};

Vipster::IO::BaseData pwi_file_parser(std::string fn, std::ifstream &file)
{
    Vipster::IO::PWData d;
    d.mol.setName(fn);
//    Vipster::IO::PWParam* p = dynamic_cast<Vipster::IO::PWParam*>(&d.data);

    typedef std::istreambuf_iterator<char> iter;
    boost::spirit::multi_pass<iter> first = boost::spirit::make_default_multi_pass(iter(file));

    pwi_parse_grammar<boost::spirit::multi_pass<iter>> grammar;

    qi::phrase_parse(first, boost::spirit::make_default_multi_pass(iter()), grammar, qi::blank, d);

//    cout << *p << endl;
    return d;
}

const Vipster::IOPlugin Vipster::IO::PWInput =
{
    "pwi",
    "pwi",
    "pwi",
    &pwi_file_parser,
    nullptr
};
