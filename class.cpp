#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <numeric>

#include <iostream>


// only made solution into a class to minimize the shit I have to pass around, just to minimize the number of columns I use
class Solution
{
public:
    Solution(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, const double targetBuySell, const std::size_t maxIter) noexcept
        : Solution(targetWeights, assetValues, targetBuySell, maxIter, std::accumulate(assetValues.begin(), assetValues.end(), 0))
    {}

    Solution(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, const double targetBuySell, const std::size_t maxIter, const double totalPortfolioValue) noexcept
        : targetWeights(targetWeights), assetValues(assetValues), targetBuySell(targetBuySell), maxIter(maxIter), totalPortfolioValue(totalPortfolioValue), isSell(targetBuySell < 0)
    {}

    Solution() = delete;
    Solution(const Solution&) = delete;
    Solution& operator=(const Solution&) = delete;
    Solution(Solution&&) = delete;
    Solution& operator=(Solution&&) = delete;

    std::vector<double> getResult(void) const noexcept;

private:
    const std::vector<double>& targetWeights;
    const std::vector<double>& assetValues;
    const double targetBuySell;
    const std::size_t maxIter;
    const double totalPortfolioValue;
    const bool isSell;

    bool currencyEquality(const double a, const double b) const noexcept;
    bool isWeightZero(const double weight) const noexcept;

    bool canDoFullSingleOpRebalance(void /* targetWeights, assetValues, totalPortfolioValue, targetBuySell */) const noexcept;
    std::size_t getMostOverWeightIndex(/* targetWeights, */ const std::vector<double>&currAssetValues, const double currPortfolioValue) const noexcept;
    std::size_t getMostUnderWeightIndex(/* targetWeights, */ const std::vector<double>&currAssetValues, const double currPortfolioValue) const noexcept;
    std::vector<double> singleOperationTransaction(void /* targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter */) const noexcept;
    std::vector<double> zeroWeightHandler(/* targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter, */ const std::size_t numberOfZeroWeights) const noexcept;
    std::vector<double> sellZeroWeightAsset(const std::vector<double>& zeroWeightAssetValues, const std::vector<std::size_t>& zeroWeightOriginalIndexNumbers, const double targetSell) const noexcept;
};


template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "[ ";
    for (std::size_t i = 0; i < vec.size() - 1; i += 1)
    {
        os << vec[i] << ", ";
    }

    os << vec[vec.size() - 1] << "]";
    return os;
}


int main(void)
{
    std::vector<double> targetWeights = {
        65.000 / 100,
        0.000 / 100,
        12.000 / 100,
        10.000 / 100,
        0.000 / 100,
        0.000 / 100,
        13.000 / 100
    };

    std::vector<double> assetValues = {
        23672.72,
        15918.04,
        6708.47,
        3659.90,
        5608.27,
        14447.1967749764,
        11317.59
    };

    const double targetBuySell = -20000;
    const auto a = Solution(targetWeights, assetValues, targetBuySell, abs(targetBuySell));

    std::cout << Solution(targetWeights, assetValues, targetBuySell, abs(targetBuySell)).getResult() << std::endl;

    return 0;
}




bool floatingPointEquality(double a, double b, double epsilon) noexcept
{
    return std::abs(a - b) <= epsilon;
}




std::vector<double> Solution::getResult(void) const noexcept
{
    return this->singleOperationTransaction();
}





bool Solution::currencyEquality(double a, double b) const noexcept
{
    // min value in a dollar is 0.01 (aka a penny). epsilon of a tenth of a penny should be good enough.
    constexpr double epsilon = 0.001;
    return floatingPointEquality(a, b, epsilon);
}



bool Solution::isWeightZero(double weight) const noexcept
{
    // I only really want weights with three decimal digits max. (ex. 5.535%). so if weight is less than or eq to 0.0001% I want to treat as zero
    constexpr double epsilon = 1e-6;
    return floatingPointEquality(weight, 0.0, epsilon);
}




/**
 * @brief the number of target weights that would return true for
 * Solution::isWeightZero SHALL be ZERO. Target weights are not
 * allowed and WILL break this function.
 * 
 */
bool Solution::canDoFullSingleOpRebalance(void /* targetWeights, assetValues, totalPortfolioValue, targetBuySell */) const noexcept
{
    const auto changeInTotalPortfolioValueNeededForAssetValueToBecomeTheTargetValue = [](
        double weight,
        double assetValue,
        double portfolioValue
    ){return (assetValue / weight) - portfolioValue;};

    if (this->isSell)
    {
        double minPortfolioSellForFullOnlySellRebalance = this->targetWeights[0];
        for (std::size_t i = 1; i < this->targetWeights.size(); i += 1)
        {
            minPortfolioSellForFullOnlySellRebalance = std::max(
                minPortfolioSellForFullOnlySellRebalance,
                changeInTotalPortfolioValueNeededForAssetValueToBecomeTheTargetValue(
                    this->targetWeights[i],
                    this->assetValues[i],
                    this->totalPortfolioValue
                )
            );
        }
        const double targetSell = targetBuySell; // is a negative value

        
    }
}





std::size_t Solution::getMostOverWeightIndex(/* targetWeights, */ const std::vector<double>&currAssetValues, const double currPortfolioValue) const noexcept
{

}





std::size_t Solution::getMostUnderWeightIndex(/* targetWeights, */ const std::vector<double>&currAssetValues, const double currPortfolioValue) const noexcept
{

}





std::vector<double> Solution::singleOperationTransaction(void /* targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter */) const noexcept
{

}





std::vector<double> Solution::zeroWeightHandler(/* targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter, */ const std::size_t numberOfZeroWeights) const noexcept
{

}





std::vector<double> Solution::sellZeroWeightAsset(const std::vector<double>& zeroWeightAssetValues, const std::vector<std::size_t>& zeroWeightOriginalIndexNumbers, const double targetSell) const noexcept
{

}




