#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <numeric>

#include <iostream>


// only made solution into a class to minimize the shit I have to pass around, just to minimize
// the number of columns I use and to do some optimization in single operation transaction
/**
 * SUM(targetWeights) must eq 100%
 * targetWeights.size() == assetValues.size()
 * targetWeights.size() >= 2 (that'd be a waste if it was less than 2...)
 * maxIter > 0
 * targetBuySell should != 0.0 (that'd be a waste if it was..,)
 * totalPortfolioValue = SUM(assetValues)
 */
class Solution
{
public:
    Solution(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, const double totalPortfolioValue, const double targetBuySell, const std::size_t maxIter) noexcept
        : targetWeights(targetWeights), assetValues(assetValues), totalPortfolioValue(totalPortfolioValue), targetBuySell(targetBuySell), maxIter(maxIter), isSell(targetBuySell < 0), assetCount(assetValues.size())
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
    const std::size_t assetCount;
    const bool isSell;

    bool canDoFullSingleOpRebalance(void /* targetWeights, assetValues, totalPortfolioValue, targetBuySell */) const noexcept;
    std::size_t getMostOverWeightIndex(/* targetWeights, */ const std::vector<double>& currAssetValues, const double currPortfolioValue) const noexcept;
    std::size_t getMostUnderWeightIndex(/* targetWeights, */ const std::vector<double>& currAssetValues, const double currPortfolioValue) const noexcept;
    std::vector<double> singleOperationTransaction(void /* targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter */) const noexcept;
    std::vector<double> zeroWeightHandler(/* targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter, */ const std::size_t numberOfZeroWeights) const noexcept;
    std::vector<double> sellZeroWeightAsset(const std::vector<double>& zeroWeightAssetValues, const std::vector<std::size_t>& zeroWeightOriginalIndexNumbers, const double targetSell) const noexcept;

    static bool isWeightZero(const double weight) noexcept;
    static bool isCurrencyEqual(const double a, const double b) noexcept;
    static double rebalanceAssetValue(const double targetWeight, const double assetValue, const double portfolioValue) noexcept;
};

bool isFloatingPointEqual(const double a, const double b, const double epsilon) noexcept;

// for testing only, will be deleted in the future
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
    const std::vector<double> targetWeights = {
        65.000 / 100,
        0.000 / 100,
        12.000 / 100,
        10.000 / 100,
        0.000 / 100,
        0.000 / 100,
        13.000 / 100
    };
    const std::vector<double> assetValues = {
        23672.72,
        15918.04,
        6708.47,
        3659.90,
        5608.27,
        14447.1967749764,
        11317.59
    };
    const double totalPortfolioValue = std::accumulate(assetValues.begin(), assetValues.end(), 0);
    const double targetBuySell = -20000;
    const std::size_t maxIter = std::abs(targetBuySell) + 1;

    std::cout << Solution(targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter).getResult() << std::endl;

    return 0;
}




bool isFloatingPointEqual(const double a, const double b, const double epsilon) noexcept
{
    return std::abs(a - b) <= epsilon;
}




double Solution::rebalanceAssetValue(const double targetWeight, const double assetValue, const double portfolioValue) noexcept
{
    const double idealValue = portfolioValue * targetWeight;
    return idealValue - assetValue;
}




bool Solution::isCurrencyEqual(double a, double b) noexcept
{
    // min value in a dollar is 0.01 (aka a penny). epsilon is a tenth of a penny. That should be good enough.
    constexpr double epsilon = 0.001;
    return isFloatingPointEqual(a, b, epsilon);
}



bool Solution::isWeightZero(double weight) noexcept
{
    // if a weight is less than 1e-9, that means in order for the portfolio to be able to do a
    // full single op rebalance would take more than 1 billion dollars... (look at canDoFullSingleOpRebalance func)
    constexpr double epsilon = 1e-9;
    return isFloatingPointEqual(weight, 0.0, epsilon);
}




std::vector<double> Solution::getResult(void) const noexcept
{
    const bool isSellWholePortfolio = -this->targetBuySell > this->totalPortfolioValue || Solution::isCurrencyEqual(-this->targetBuySell, this->totalPortfolioValue);
    if (isSellWholePortfolio)
    {
        std::vector<double> fullSell(this->assetCount);
        std::transform(assetValues.begin(), assetValues.end(), fullSell.begin(), [](double value){return -value;});
        return fullSell;
    }

    const std::size_t numberOfZeroWeights = std::count_if(this->targetWeights.begin(), this->targetWeights.end(), Solution::isWeightZero);
    if (numberOfZeroWeights > 0)
        return this->zeroWeightHandler(numberOfZeroWeights);

    return this->singleOperationTransaction();
}




/**
 * @brief the number of target weights that would return true for
 * Solution::isWeightZero SHALL be ZERO. Target weights are not
 * allowed and WILL break this function.
 * 
 */
bool Solution::canDoFullSingleOpRebalance(void /* targetWeights, assetValues, totalPortfolioValue, targetBuySell */) const noexcept
{
    // for some reason when I use the ternary operator, std::min and std::max don't know which function to overload; use static_cast to dictate this.
    using extremum_t = const double& (*)(const double&, const double&);
    const extremum_t funcExtremum = (this->isSell) ? static_cast<extremum_t>(std::min<double>) : static_cast<extremum_t>(std::max<double>);

    const double minChangeForFullSingleOpRebalance = ([&]()
    {
        // if asset is over weight this value will be positive and if underweight this will be negative
        // balance means have target weight or  target value, which is = total portfolio value * target weight
        const auto changeInPortfolioValueForCurrAssetToBeBalanced = [&](std::size_t i)
            {return (this->assetValues[i] / this->targetWeights[i]) - this->totalPortfolioValue;};

        double result = changeInPortfolioValueForCurrAssetToBeBalanced(0);
        for (std::size_t i = 1; i < this->assetCount; i += 1)
        {
            result = funcExtremum(result, changeInPortfolioValueForCurrAssetToBeBalanced(i));
        }

        return result;
    })();

    if (this->isSell)
    {
        const double targetSell = this->targetBuySell; // is a negative value
        return -targetSell > -minChangeForFullSingleOpRebalance || Solution::isCurrencyEqual(targetSell, minChangeForFullSingleOpRebalance);
    }

    const double targetBuy = this->targetBuySell;
    return targetBuy > minChangeForFullSingleOpRebalance || Solution::isCurrencyEqual(targetBuy, minChangeForFullSingleOpRebalance);
}





std::size_t Solution::getMostOverWeightIndex(/* targetWeights, */ const std::vector<double>& currAssetValues, const double currPortfolioValue) const noexcept
{
    std::size_t overWeightIndex = 0;
    double overWeightValue = Solution::rebalanceAssetValue(this->targetWeights[0], currAssetValues[0], currPortfolioValue);
    for (std::size_t i = 1; i < this->assetCount; i += 1)
    {
        const double iValue = Solution::rebalanceAssetValue(this->targetWeights[i], currAssetValues[i], currPortfolioValue);
        // looking for the asset with the most negative rebalance change
        if (iValue < overWeightValue)
        {
            overWeightIndex = i;
            overWeightValue = iValue;
        }
    }

    return overWeightIndex;
}





std::size_t Solution::getMostUnderWeightIndex(/* targetWeights, */ const std::vector<double>& currAssetValues, const double currPortfolioValue) const noexcept
{
    std::size_t underWeightIndex = 0;
    double underWeightValue = Solution::rebalanceAssetValue(this->targetWeights[0], currAssetValues[0], currPortfolioValue);
    for (std::size_t i = 1; i < this->assetCount; i += 1)
    {
        const double iValue = Solution::rebalanceAssetValue(this->targetWeights[i], currAssetValues[i], currPortfolioValue);
        // looking for the asset with the most positive rebalance change
        if (iValue > underWeightValue)
        {
            underWeightIndex = i;
            underWeightValue = iValue;
        }
    }

    return underWeightIndex;
}





std::vector<double> Solution::singleOperationTransaction(void /* targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter */) const noexcept
{

}





std::vector<double> Solution::zeroWeightHandler(const std::size_t numberOfZeroWeights /*, targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter */) const noexcept
{
}





std::vector<double> Solution::sellZeroWeightAsset(const std::vector<double>& zeroWeightAssetValues, const std::vector<std::size_t>& zeroWeightOriginalIndexNumbers, const double targetSell) const noexcept
{

}




