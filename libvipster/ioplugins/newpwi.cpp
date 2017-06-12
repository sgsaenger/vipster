#include "ioplugins/pwinput.h"
#define BOOST_SPIRIT_DEBUG
#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/phoenix.hpp>
//#include <boost/fusion/include/adapt_struct_named.hpp>
//#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/fusion/include/std_pair.hpp>

namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;
using namespace std;
using namespace Vipster;
using namespace qi;
using pwmap = std::map<std::string,std::string>;

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::IO::PWData,
        data
//        mol
)

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
        : grammar<Iterator, IO::PWData(), blank_type>
{
    pwi_parse_grammar(): pwi_parse_grammar::base_type(file)
    {
        file = param > mol;
        file.name("PWScf input file");
        /* Namelists:
         * &NAME
         * key=value, key=value
         * key=value
         * !comment
         * /
         */
        nl_names.add
              ("control", &IO::PWParam::control)
              ("system", &IO::PWParam::system)
              ("electrons", &IO::PWParam::electrons)
              ("ions", &IO::PWParam::ions)
              ("cell", &IO::PWParam::cell);
        param = ('&' > no_case[nl_names[_a = _1]] > eol >
                 nl[phx::bind([](pwmap IO::PWParam::* v, IO::PWParam& p, pwmap m){p.*v = m;cout << "AHA!\n";}, _a, _val, _1)]
                 > eol > '/')
                % *eol;
        nl = nl_entry % (eol|',');
        nl_entry = key > qi::lit('=') > value;
        value = fstr | key;
        key = +(qi::graph - qi::char_("=/,'"));
        fstr = qi::char_('\'') >> qi::lexeme[*(qi::char_ - '\'')] > qi::char_('\'');
        /* Cards:
         * NAME {options}
         * #comment
         * !comment
         * <content>
         */
        mol = *(no_case[card_names[_a = _1]] > lazy(*_a)(_val));
        card_names.add
                ("atomic_positions",nullptr)
                ("atomic_species",nullptr)
                ("cell_parameters",nullptr)
                ("k_points",nullptr);
        on_error<fail>(
                file,
                cout
                << phx::val("Error! Expecting ")
                << _4
                << phx::val(" here: \"")
                << phx::construct<std::string>(_3, _2)
                << phx::val("\"")
                << endl
        );
        BOOST_SPIRIT_DEBUG_NODE(key);
        BOOST_SPIRIT_DEBUG_NODE(value);
        BOOST_SPIRIT_DEBUG_NODE(fstr);
    }
    rule<Iterator, IO::PWData(), blank_type> file;
    //PWParam
    rule<Iterator, IO::PWParam(), locals<pwmap IO::PWParam::*>, blank_type> param;
    rule<Iterator, pwmap(), blank_type> nl;
    symbols<char, pwmap IO::PWParam::*> nl_names;
    rule<Iterator, pair<std::string,std::string>(), blank_type> nl_entry;
    rule<Iterator, std::string()> key;
    rule<Iterator, std::string()> value;
    rule<Iterator, std::string()> fstr;
    //Molecule
    rule<Iterator, Molecule(), locals<rule<Iterator, blank_type>*>, blank_type> mol;
    symbols<char, rule<Iterator, blank_type>*> card_names;
};

Vipster::IO::BaseData pwi_file_parser(std::string fn, std::ifstream &file)
{
    // set up data
    Vipster::IO::PWData d;
    d.mol.setName(fn);

    // set up iterators
    typedef std::istreambuf_iterator<char> iter;
    boost::spirit::multi_pass<iter> first = boost::spirit::make_default_multi_pass(iter(file));
    boost::spirit::multi_pass<iter> last = boost::spirit::make_default_multi_pass(iter());

    pwi_parse_grammar<boost::spirit::multi_pass<iter>> grammar;

    phrase_parse(first, last, grammar, blank, d);

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
