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
        Vipster::IO::PWParam,
        control,
        system,
        electrons,
        ions,
        cell
)

template<typename Iterator>
struct pwi_parse_grammar
        : qi::grammar<Iterator, shared_ptr<Vipster::IO::PWParam>(), qi::blank_type>
{
    pwi_parse_grammar(): pwi_parse_grammar::base_type(file)
    {
        key = *(qi::graph - qi::char_("=/,'"));
        fstr = qi::char_('\'') >> qi::lexeme[*(qi::char_ - '\'')] > qi::char_('\'');
        value = fstr | key;
        nl_entry = key >> qi::lit('=') >> value;
        namelist = *qi::eol
                   >> nl_entry % (qi::eol | ',')
                   >> *qi::eol >> '/' >> *qi::eol;
        param = qi::no_case["&control"] > namelist >
                qi::no_case["&system"] > namelist >
                qi::no_case["&electrons"] > namelist >
                -(qi::no_case["&ions"] > namelist) >
                -(qi::no_case["&cell"] > namelist);
        file = param;
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
        BOOST_SPIRIT_DEBUG_NODE(file);
        BOOST_SPIRIT_DEBUG_NODE(param);
        BOOST_SPIRIT_DEBUG_NODE(namelist);
        BOOST_SPIRIT_DEBUG_NODE(nl_entry);
        BOOST_SPIRIT_DEBUG_NODE(key);
        BOOST_SPIRIT_DEBUG_NODE(fstr);
        BOOST_SPIRIT_DEBUG_NODE(value);
    }
//    qi::rule<Iterator, Vipster::IOData(), qi::blank_type> file;
    qi::rule<Iterator, shared_ptr<Vipster::IO::PWParam>(), qi::blank_type> file;
    qi::rule<Iterator, Vipster::IO::PWParam(), qi::blank_type> param;
    qi::rule<Iterator, map<string, string>(), qi::blank_type> namelist;
    qi::rule<Iterator, pair<string, string>(), qi::blank_type> nl_entry;
    qi::rule<Iterator, string(), qi::blank_type> key;
    qi::rule<Iterator, string(), qi::blank_type> value;
    qi::rule<Iterator, string(), qi::blank_type> fstr;
    qi::rule<Iterator, Vipster::Molecule(), qi::blank_type> mol;
};

Vipster::IOData pwi_file_parser(std::string fn, std::ifstream &file)
{
    Vipster::IOData d{
        Vipster::Molecule{fn,1},
        shared_ptr<Vipster::IO::PWParam>(new Vipster::IO::PWParam)
    };
    Vipster::IO::PWParam* p = dynamic_cast<Vipster::IO::PWParam*>(d.data.get());

    typedef std::istreambuf_iterator<char> iter;
    boost::spirit::multi_pass<iter> first = boost::spirit::make_default_multi_pass(iter(file));

    pwi_parse_grammar<boost::spirit::multi_pass<iter>> grammar;

    qi::phrase_parse(first, boost::spirit::make_default_multi_pass(iter()), grammar, qi::blank, d.data);

    cout << *p << endl;
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
