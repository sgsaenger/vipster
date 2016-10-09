struct Atom{
    std::string name;
    Vec coord;
    float charge;
    std::array<bool,3> fix;
    bool hidden;
};
