#include "HelperRepos/jsArray.h"
#include <limits>
#include <cmath>
#include <numeric>


#include <iostream>


// only made solution into a class to minimize the shit I have to pass around, just to minimize
// the number of columns I use and to do some optimization in single operation transaction

/**
 * SUM(targetWeights) must eq 100%
 * targetWeights.size() == assetValues.size()
 * totalPortfolioValue = SUM(assetValues)
 * maxIter > 0
 */
class Solution
{
public:
    Solution(const JSArray<double>& targetWeights, const JSArray<double>& assetValues, const double targetBuySell, const std::size_t maxIter) noexcept
        : Solution(
            targetWeights,
            assetValues,
            std::accumulate(assetValues.begin(), assetValues.end(), 0),
            std::count_if(targetWeights.begin(), targetWeights.end(), Solution::isWeightZero),
            targetBuySell,
            maxIter
        )
    {}

    JSArray<double> getResult(void) const noexcept;


    Solution() = delete;
    Solution(const Solution&) = delete;
    Solution& operator=(const Solution&) = delete;
    Solution(Solution&&) = delete;
    Solution& operator=(Solution&&) = delete;

private:
    Solution(const JSArray<double>& targetWeights, const JSArray<double>& assetValues, const double totalPortfolioValue, const std::size_t numberOfZeroWeights, const double targetBuySell, const std::size_t maxIter) noexcept
        : targetWeights(targetWeights),
        assetValues(assetValues),
        totalPortfolioValue(totalPortfolioValue),
        numberOfZeroWeights(numberOfZeroWeights),
        targetBuySell(targetBuySell),
        maxIter(maxIter),
        assetCount(assetValues.size()),
        isSell(targetBuySell < 0)
    {}

    const JSArray<double>& targetWeights;
    const JSArray<double>& assetValues;
    const double totalPortfolioValue;
    const std::size_t numberOfZeroWeights;
    const double targetBuySell;
    const std::size_t maxIter;
    const std::size_t assetCount;
    const bool isSell;

    bool canDoFullSingleOpRebalance(void /* targetWeights, assetValues, totalPortfolioValue, targetBuySell */) const noexcept;
    std::size_t getMostOverWeightIndex(/* targetWeights, */ const JSArray<double>& currAssetValues, const double currPortfolioValue) const noexcept;
    std::size_t getMostUnderWeightIndex(/* targetWeights, */ const JSArray<double>& currAssetValues, const double currPortfolioValue) const noexcept;
    JSArray<double> singleOperationTransaction(void /* targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter */) const noexcept;
    JSArray<double> zeroWeightHandler(/* ,targetWeights, assetValues, totalPortfolioValue, numberOfZeroWeights, targetBuySell, maxIter */) const noexcept;
    JSArray<double> sellZeroWeightAsset(const JSArray<double>& zeroWeightAssetValues, const JSArray<std::size_t>& zeroWeightOriginalIndexNumbers, const double targetSell) const noexcept;

    static bool isWeightZero(const double weight) noexcept;
    static bool isCurrencyEqual(const double a, const double b) noexcept;
    static double rebalanceAssetValue(const double targetWeight, const double assetValue, const double portfolioValue) noexcept;
    static double getAssetValueChangeForEqRebalance(const double weight0, const double value0, const double weight1, const double value1, const double portfolioValue) noexcept;
};

bool isFloatingPointEqual(const double a, const double b, const double epsilon) noexcept;

// for testing only, will be deleted in the future
template<typename T>
std::ostream& operator<<(std::ostream& os, const JSArray<T>& vec) {
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
    const JSArray<double> targetWeights = {
        65.000 / 100,
        5.000 / 100,
        6.000 / 100,
        5.000 / 100,
        5.000 / 100,
        6.000 / 100,
        8.000 / 100
    };
    const JSArray<double> assetValues = {
        23672.72,
        15918.04,
        6708.47,
        3659.90,
        5608.27,
        14447.1967749764,
        11317.59
    };
    const double targetBuySell = 71961.4 + 3700;

    auto calculator = Solution(targetWeights, assetValues, targetBuySell, 100);
    std::cout << calculator.getResult() << '\n';


    return 0;
}




bool isFloatingPointEqual(const double a, const double b, const double epsilon) noexcept
{
    return std::abs(a - b) <= epsilon;
}





bool Solution::isWeightZero(double weight) noexcept
{
    // if a weight is less than 1e-9, that means in order for the portfolio to be able to do a
    // full single op rebalance would take more than 1 billion dollars... (look at canDoFullSingleOpRebalance func)
    constexpr double epsilon = 1e-9;
    return isFloatingPointEqual(weight, 0.0, epsilon);
}




bool Solution::isCurrencyEqual(double a, double b) noexcept
{
    // min value in a dollar is 0.01 (aka a penny). epsilon is a tenth of a penny. That should be good enough.
    constexpr double epsilon = 0.001;
    return isFloatingPointEqual(a, b, epsilon);
}



double Solution::rebalanceAssetValue(const double targetWeight, const double assetValue, const double portfolioValue) noexcept
{
    const double idealValue = portfolioValue * targetWeight;
    return idealValue - assetValue;
}



double Solution::getAssetValueChangeForEqRebalance(const double weight0, const double value0, const double weight1, const double value1, const double portfolioValue) noexcept
{
    /**
     * Calculation Notes:
     * rebalance(pv, tw, av) = pv * tw - av;
     * rebalance(pv + x, tw0, av0 + x) = rebalance(pv + x, tw1, av1)  // where x is the change in av0 so values are equal. (if av0 is most underweight x will always be > 0, if av0 is most overweight x will always be < 0)
     * (pv + x) * tw0 - (av0 + x) = (pv + x) * tw1 - av1
     * pv*tw0 + x*tw0 - av0 - x = pv*tw1 + x*tw1 - av1
     * x*(tw0 - tw1 - 1) + pv*tw0 - av0 = pv*tw1 - av1
     * x*(tw0 - tw1 - 1) = pv*(tw1 - tw0) + av0 - av1
     * x = (pv*(tw1 - tw0) + av0 - av1) / (tw0 - tw1 - 1)
     * 
     */
    return (portfolioValue * (weight1 - weight0) + value0 - value1) / (weight0 - weight1 - 1);
}




JSArray<double> Solution::getResult(void) const noexcept
{
    if (Solution::isCurrencyEqual(this->targetBuySell, 0.0))
        return JSArray<double>(this->assetCount, 0.0);

    if (this->assetCount <= 1)
        return JSArray<double>(this->assetCount, this->targetBuySell);

    const bool isSellWholePortfolio = -this->targetBuySell > this->totalPortfolioValue || Solution::isCurrencyEqual(-this->targetBuySell, this->totalPortfolioValue);
    if (isSellWholePortfolio)
        return this->assetValues.map([](double val){return -val;});

    if (this->numberOfZeroWeights > 0)
        return this->zeroWeightHandler();

    if (this->canDoFullSingleOpRebalance())
        return this->targetWeights.map(
            [&](double weight, std::size_t i)
            {
                return Solution::rebalanceAssetValue(weight, this->assetValues[i], this->totalPortfolioValue);
            }
        );

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
    const double minChangeForFullSingleOpRebalance = ([&]()
    {
        // for some reason when I use the ternary operator, std::min and std::max don't know which function to overload; use static_cast to dictate this.
        using extremum_t = const double& (*)(const double&, const double&);
        const extremum_t getChangeNeededForFullRebalance = (this->isSell) ? static_cast<extremum_t>(std::min<double>) : static_cast<extremum_t>(std::max<double>);

        // if asset is over weight this value will be positive and if underweight this will be negative
        // balance means have target weight or  target value, which is = total portfolio value * target weight
        const auto changeInPortfolioValueForCurrAssetToBeBalanced = [&](std::size_t i)
            {return (this->assetValues[i] / this->targetWeights[i]) - this->totalPortfolioValue;};

        double result = changeInPortfolioValueForCurrAssetToBeBalanced(0);
        for (std::size_t i = 1; i < this->assetCount; i += 1)
        {
            result = getChangeNeededForFullRebalance(result, changeInPortfolioValueForCurrAssetToBeBalanced(i));
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





std::size_t Solution::getMostOverWeightIndex(/* targetWeights, */ const JSArray<double>& currAssetValues, const double currPortfolioValue) const noexcept
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





std::size_t Solution::getMostUnderWeightIndex(/* targetWeights, */ const JSArray<double>& currAssetValues, const double currPortfolioValue) const noexcept
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





JSArray<double> Solution::singleOperationTransaction(void /* targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter */) const noexcept
{
    const double changePacket = this->targetBuySell / this->maxIter;
    const auto assetSelector = (this->isSell) ? &Solution::getMostOverWeightIndex : &Solution::getMostUnderWeightIndex;
    const std::size_t initTargetAssetIndex = (this->*assetSelector)(assetValues, totalPortfolioValue);
    const double initChange = ([&]() -> double
    {
        using extremum_t = const double& (*)(const double&, const double&);
        const extremum_t getSmallestAbsChangeForEqRebalance = (this->isSell) ? static_cast<extremum_t>(std::max<double>) : static_cast<extremum_t>(std::min<double>);
        // want the smallest change to keep the initTargetAssetIndex still tied for the most over/under weight asset

        const double weight0 = targetWeights[initTargetAssetIndex];
        const double value0 = assetValues[initTargetAssetIndex];
        const auto getValue = [&](std::size_t i) -> double
        {
            const double weight1 = this->targetWeights[i];
            const double value1 = this->assetValues[i];
            return Solution::getAssetValueChangeForEqRebalance(weight0, value0, weight1, value1, this->totalPortfolioValue);
        };

        double result = this->targetBuySell;
        for (std::size_t i = 0; i < initTargetAssetIndex; i += 1)
        {
            result = getSmallestAbsChangeForEqRebalance(result, getValue(i));
        }

        for (std::size_t i = initTargetAssetIndex + 1; i < this->assetCount; i += 1)
        {
            result = getSmallestAbsChangeForEqRebalance(result, getValue(i));
        }

        return result;
    })();
    JSArray<double> currAssetValues = this->assetValues;
    double currPortfolioValue = totalPortfolioValue + initChange;
    const double targetPortfolioValue = this->totalPortfolioValue + this->targetBuySell;
    const std::size_t remainingIter = (targetBuySell - initChange) / changePacket;

    currAssetValues[initTargetAssetIndex] += initChange;
    for (std::size_t i = 0; i < remainingIter; i += 1, currPortfolioValue += changePacket)
    {
        currAssetValues[(this->*assetSelector)(currAssetValues, currPortfolioValue)] += changePacket;
    }
    // currPortfolioValue will always be <= targetPortfolioValue.
    currAssetValues[(this->*assetSelector)(currAssetValues, currPortfolioValue)] += targetPortfolioValue - currPortfolioValue;

    return ([&]() -> JSArray<double>
    {
        JSArray<double>& diff = currAssetValues;
        for (std::size_t i = 0; i < diff.size(); i += 1)
        {
            diff[i] -= this->assetValues[i];
        }

        return diff;
    })();
}





JSArray<double> Solution::zeroWeightHandler(/* ,targetWeights, assetValues, totalPortfolioValue, numberOfZeroWeights, targetBuySell, maxIter */) const noexcept
{
    const auto alignResultsToOriginalIndex = [](
        const JSArray<double>& subsetValues,
        const JSArray<std::size_t>& originalSubsetElementIndex,
        std::size_t superSetSize
    ) -> JSArray<double>
    {
        JSArray<double> aligned(superSetSize);
        for (std::size_t i = 0, k = 0; i < superSetSize; i += 1)
        {
            aligned[i] = ([&]() -> double
            {
                if (k < originalSubsetElementIndex.size())
                {
                    const std::size_t originalIndex = originalSubsetElementIndex[k];
                    if (originalIndex == i)
                    {
                        const double result = subsetValues[k];
                        k += 1;
                        return result;
                    }
                }

                return 0;
            })();
        }

        return aligned;
    };

    const JSArray<std::size_t> indexNumbersNonZeroWeight = ([&]()
    {
        JSArray<std::size_t> result(this->assetCount - this->numberOfZeroWeights);
        for (std::size_t i = 0, k = 0; i < this->assetCount; i += 1)
        {
            const double weight = this->targetWeights[i];
            if (!Solution::isWeightZero(weight))
            {
                result[k] = i;
                k += 1;
            }
        }

        return result;
    })();
    const JSArray<double> targetNonZeroWeights = indexNumbersNonZeroWeight
        .map([&](std::size_t valIndex){return this->targetWeights[valIndex];});
    const JSArray<double> assetValuesWithTargetNonZeroWeights = indexNumbersNonZeroWeight
        .map([&](std::size_t valIndex){return this->assetValues[valIndex];});
    const double sumOfNonZeroWeightAssets = assetValuesWithTargetNonZeroWeights
        .reduce([](double sum, double val){return sum + val;}, 0.0);

    if (!this->isSell)
    {
        return alignResultsToOriginalIndex(
            Solution(
                targetNonZeroWeights,
                assetValuesWithTargetNonZeroWeights,
                sumOfNonZeroWeightAssets,
                0,
                this->targetBuySell,
                this->maxIter
            ).getResult(),
            indexNumbersNonZeroWeight,
            this->assetCount
        );
    }

    const double targetSell = targetBuySell;
    const std::size_t numberOfZeroWeights = this->assetCount - indexNumbersNonZeroWeight.size();
    const JSArray<std::size_t> indexNumbersZeroWeight = ([&]()
    {
        JSArray<std::size_t> result(this->numberOfZeroWeights);
        for (std::size_t i = 0, k = 0; i < this->assetCount; i += 1)
        {
            const double weight = this->targetWeights[i];
            if (Solution::isWeightZero(weight))
            {
                result[k] = i;
                k += 1;
            }
        }

        return result;
    })();
    const double sumOfZeroWeightAssetValues = this->totalPortfolioValue - sumOfNonZeroWeightAssets;
    const bool canSellMoreThanJustAllTheZeroWeightAssets = -targetSell > sumOfZeroWeightAssetValues;
    const bool canOnlySellSomeZeroWeightAssets = sumOfZeroWeightAssetValues > -targetSell;

    if (canSellMoreThanJustAllTheZeroWeightAssets)
    {
        const std::size_t newMaxIter = ([&]() -> std::size_t
        {
            const std::size_t originalMaxIter = this->maxIter;
            const std::size_t originalChangePacket = totalPortfolioValue / originalMaxIter;
            const std::size_t equivalentItersDoneBySellingAllZeroWeightAssets = sumOfZeroWeightAssetValues / originalChangePacket;
            return originalMaxIter - equivalentItersDoneBySellingAllZeroWeightAssets;
        })();
        const JSArray<double> changesToNonZeroAssetValues = Solution(
                targetNonZeroWeights,
                assetValuesWithTargetNonZeroWeights,
                sumOfNonZeroWeightAssets,
                0,
                this->targetBuySell + sumOfZeroWeightAssetValues,
                newMaxIter
            ).getResult();

        return ([&]() -> JSArray<double>
        {
            JSArray<double> result = alignResultsToOriginalIndex(
                changesToNonZeroAssetValues,
                indexNumbersNonZeroWeight,
                targetWeights.size()
            );
            for (const std::size_t zeroIndex : indexNumbersZeroWeight)
            {
                result[zeroIndex] = -assetValues[zeroIndex];
            }

            return result;
        })();
    }
}





JSArray<double> Solution::sellZeroWeightAsset(const JSArray<double>& zeroWeightAssetValues, const JSArray<std::size_t>& zeroWeightOriginalIndexNumbers, const double targetSell) const noexcept
{
    // I need to be able to sort and still know what the element's
    // original index was in order to put them back into their
    // original order.
    struct ValueIndexPair
    {
        ValueIndexPair(double value, std::size_t index)
            : value(value), index(index)
        {}

        double value;
        int index;
    };
    const double targetSell = this->targetBuySell;

    if (this->numberOfZeroWeights == 1) return {zeroWeightAssetValues[0] + targetSell};

    JSArray<ValueIndexPair> zeroWeightAssetsSortedValuesDescending = zeroWeightAssetValues
        .map([](double value, std::size_t index){return ValueIndexPair(value, index);})
        .sort([](const auto& a, const auto& b){return a.value > b.value});
    std::size_t groupWithSameLargestValueSize = ([&]() -> std::size_t
    {
        std::size_t groupSize = 1;
        const double largestValue = zeroWeightAssetsSortedValuesDescending[0].value;
        for (;
            groupSize < zeroWeightAssetsSortedValuesDescending.size() &&
                Solution::isCurrencyEqual(largestValue, zeroWeightAssetsSortedValuesDescending[groupSize].value);
            groupSize += 1
        ){}

        return groupSize;
    })();
}




