#include "ioplugins/pwinput.h"
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/adapt_adt.hpp>

namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;

template<typename Iterator>
struct pwi_parse_grammar
        : qi::grammar<Iterator, Vipster::Molecule(), qi::blank_type>
{
    pwi_parse_grammar(): pwi_parse_grammar::base_type(mol)
    {
//        file = param >> mol;
    }
    qi::rule<Iterator, Vipster::IOData(), qi::blank_type> file;
    qi::rule<Iterator, Vipster::IO::PWParam(), qi::blank_type> param;
    qi::rule<Iterator, std::map<std::string, std::string>, qi::blank_type> namelist;
    qi::rule<Iterator, Vipster::Molecule(), qi::blank_type> mol;
};

Vipster::IOData pwi_file_parser(std::string fn, std::ifstream &file)
{
    Vipster::IOData d{
        std::make_shared<Vipster::Molecule>(fn,0),
        Vipster::IOType::None,
        std::shared_ptr<Vipster::IOBase>()
    };

    typedef std::istreambuf_iterator<char> iter;
    boost::spirit::multi_pass<iter> first = boost::spirit::make_default_multi_pass(iter(file));

    pwi_parse_grammar<boost::spirit::multi_pass<iter>> grammar;

    qi::phrase_parse(first, boost::spirit::make_default_multi_pass(iter()), grammar, qi::blank, *d.mol);

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
