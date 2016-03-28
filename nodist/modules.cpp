
// NYI ... 
//
// experiment with modules.
class Moments {
    public:
        // empty
        Moments(int order) : order(order), na_rm(false) {
            moments = Rcpp::NumericVector(1+order);
        }
        Moments(int order, bool na_rm) : order(order), na_rm(na_rm) {
            moments = Rcpp::NumericVector(1+order);
        }
        // 'pure'
        Moments(SEXP input, int order) : order(order), na_rm(false) {
            moments = wrapMoments(input, order, na_rm);
        }
        Moments(SEXP input, int order, bool na_rm) : order(order), na_rm(na_rm) {
            moments = wrapMoments(input, order, na_rm);
        }
        // append operation
        void join(const Moments& rhs) {
            moments = combineMoments(moments,rhs.moments);
        }
        // access the normalized moments
        NumericVector cent_moments(int used_df=1) {
            NumericVector vret = NumericVector(1+order);
            vret[order] = moments[0];
            vret[order-1] = moments[1];
            for (int mmm=2;mmm <= order;mmm++) {
                vret[order-mmm] = moments[mmm] / (moments[0] - used_df);
            }
            return vret;
        }
        const int order;
        const bool na_rm;
    private:
        NumericVector moments;
};

RCPP_MODULE(moment_module) {
    class_<Moments>( "Moments" )

    .constructor<int>()
    .constructor<int,bool>()
    .constructor<NumericVector,int>()
    .constructor<NumericVector,int,bool>()

    .field_readonly("order", &Moments::order)
    .field_readonly("na_rm", &Moments::na_rm)

    .method("cent_moments", &Moments::cent_moments)
    .method("%:%", &Moments::join)
    ;
}
