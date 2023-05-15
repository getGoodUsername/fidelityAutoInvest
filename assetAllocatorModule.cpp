#include "HelperRepos/jsArray.h"
#include <cmath>
#include <numeric>
#include <limits>
#include <iostream>


namespace Compute
{
    bool isFloatingPointEqual(const double a, const double b, const double epsilon) noexcept;
    bool isWeightZero(const double weight) noexcept;
    bool isCurrencyEqual(const double a, const double b) noexcept;
    double rebalanceAssetValue(const double targetWeight, const double assetValue, const double portfolioValue) noexcept;
    double getAssetValueChangeForEqRebalance(const double weight0, const double value0, const double weight1, const double value1, const double portfolioValue) noexcept;
    double getChangeInPortfolioValueForAssetValueToBecomeIdealValue(const double targetWeight, const double assetValue, const double portfolioValue) noexcept;
    double getMinTargetBuyForFullPortfolioOnlyBuyRebalance(const JSArray<double>& targetWeights, const JSArray<double>& assetValues, const double portfolioValue) noexcept;
    double getMinTargetSellForFullPortfolioOnlySellRebalance(const JSArray<double>& targetWeights, const JSArray<double>& assetValues, const double portfolioValue) noexcept;
}

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
    Solution(const JSArray<double>& targetWeights, const JSArray<double>& assetValues, const double targetBuySell) noexcept
        : Solution(
            targetWeights,
            assetValues,
            targetBuySell,
            std::abs(targetBuySell) + Compute::isCurrencyEqual(targetBuySell, 0.0)
        )
    {}

    Solution(const JSArray<double>& targetWeights, const JSArray<double>& assetValues, const double targetBuySell, const std::size_t maxIter) noexcept
        : Solution(
            targetWeights,
            assetValues,
            std::accumulate(assetValues.begin(), assetValues.end(), 0.0),
            std::count_if(targetWeights.begin(), targetWeights.end(), Compute::isWeightZero),
            targetBuySell,
            maxIter
        )
    {}

    const double totalPortfolioValue;

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
    JSArray<double> sellZeroWeightAssets(const JSArray<std::size_t>& indexNumbersZeroWeight) const noexcept;
};


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

int main(const int argc, char* argv)
{
    const JSArray<double> targetWeights = {
        65.000 / 100,
        6.000 / 100,
        0.000 / 100,
        0.000 / 100,
        0.000 / 100,
        21.000 / 100,
        9.000 / 100
    };
    const JSArray<double> assetValues = {
        23672.72,
        15918.04,
        2222,
        12222,
        24222,
        144047.1967749764,
        11317.59
    };
    const double targetBuySell = -32760;

    auto calculator = Solution(targetWeights, assetValues, targetBuySell);
    std::cout << Compute::getMinTargetBuyForFullPortfolioOnlyBuyRebalance(targetWeights, assetValues, calculator.totalPortfolioValue) << '\n';
    std::cout << Compute::getMinTargetSellForFullPortfolioOnlySellRebalance(targetWeights, assetValues, calculator.totalPortfolioValue) << '\n';
    std::cout << calculator.getResult() << '\n';

    return 0;
}



namespace Compute
{
bool isFloatingPointEqual(const double a, const double b, const double epsilon) noexcept
{
    return std::abs(a - b) <= epsilon;
}





bool isWeightZero(double weight) noexcept
{
    // if a weight is less than 1e-9, that means in order for the portfolio to be able to do a
    // full single op rebalance would take more than 1 billion dollars... (look at canDoFullSingleOpRebalance func)
    constexpr double epsilon = 1e-9;
    return Compute::isFloatingPointEqual(weight, 0.0, epsilon);
}




bool isCurrencyEqual(double a, double b) noexcept
{
    // min value in a dollar is 0.01 (aka a penny). epsilon is a tenth of a penny. That should be good enough.
    constexpr double epsilon = 0.001;
    return Compute::isFloatingPointEqual(a, b, epsilon);
}



double rebalanceAssetValue(const double targetWeight, const double assetValue, const double portfolioValue) noexcept
{
    const double idealValue = portfolioValue * targetWeight;
    return idealValue - assetValue;
}



double getAssetValueChangeForEqRebalance(const double weight0, const double value0, const double weight1, const double value1, const double portfolioValue) noexcept
{
    /**
     * Calculation Notes:
     * def rebalance(pv, tw, av) {return pv * tw - av;}
     * rebalance(pv + x, tw0, av0 + x) = rebalance(pv + x, tw1, av1)  // where x is the change in av0 so values are equal. (if av0 is most underweight x will always be > 0, if av0 is most overweight x will always be < 0)
     * (pv + x) * tw0 - (av0 + x) = (pv + x) * tw1 - av1
     * pv*tw0 + x*tw0 - av0 - x = pv*tw1 + x*tw1 - av1
     * x*(tw0 - tw1 - 1) + pv*tw0 - av0 = pv*tw1 - av1
     * x*(tw0 - tw1 - 1) = pv*(tw1 - tw0) + av0 - av1
     * x = (pv*(tw1 - tw0) + av0 - av1) / (tw0 - tw1 - 1)
     * 
     * As to why this is useful; this lets us skip a lot
     * of buy/dumb work in singleOperationTransaction. Often
     * times the most over/underweight asset's "rebalance
     * value" is much larger in magnitude the second biggest,
     * and a lot of computational cycles will be wasted trying to
     * figure out which is the next most over/under weight index
     * when it will be the same index for quite a while. We use
     * this to skip right to the point where assets are
     * swapping the position of most over/underweight asset
     * often by the just the changing by "changePacket".
     * As to why we don't use this method for every cycle,
     * when we actually get to the point where assets
     * are trading for spots for most over/under weight
     * the difference in rebalance values is not big enough
     * to warrant the extra complexity in code, but especially
     * the doubling of compute cycles. Its only worth it in the
     * beginning.
     */
    return (portfolioValue * (weight1 - weight0) + value0 - value1) / (weight0 - weight1 - 1);
}




double getChangeInPortfolioValueForAssetValueToBecomeIdealValue(const double targetWeight, const double assetValue, const double portfolioValue) noexcept
{
    /**
     * Calculation Notes:
     * av = (pv + x) * tw // where x is the change in portfolio value such that the new portfolio value (pv + x) * the target weight is av
     * av / tw = pv + x
     * x = (av / tw) - pv
     * 
     * As to why this is useful; the resulting value tells
     * us the state as to when changing the portfolio value
     * any more in the same direction (remember x can be positive
     * or negative) will result in having to change the assets
     * value too in order to keep it at its target value. If
     * we find the biggest result "x" from the set of value weight
     * pairs, this will tell us the point when enacting
     * a rebalance on all assets will at worst result in no change
     * for one asset. It therefore identifies the point as to which
     * a rebalance will result in the same operation (buy/sell) for
     * all the assets.
     */

    constexpr double positiveInfinity = std::numeric_limits<double>::infinity();
    if (Compute::isWeightZero(targetWeight))
        return (Compute::isCurrencyEqual(assetValue, 0.0)) ? 0.0 : positiveInfinity;

    return (assetValue / targetWeight) - portfolioValue;
}



double getMinTargetBuyForFullPortfolioOnlyBuyRebalance(const JSArray<double>& targetWeights, const JSArray<double>& assetValues, const double portfolioValue) noexcept
{
    return targetWeights
        .reduce(
            [&](double result, double weight, std::size_t i)
            {
                return std::max(
                    result,
                    Compute::getChangeInPortfolioValueForAssetValueToBecomeIdealValue(weight, assetValues[i], portfolioValue)
                );
            },
            0.0
        );
}



double getMinTargetSellForFullPortfolioOnlySellRebalance(const JSArray<double>& targetWeights, const JSArray<double>& assetValues, const double portfolioValue) noexcept
{
    return targetWeights
        .reduce(
            [&](double result, double weight, std::size_t i)
            {
                return std::min(
                    result,
                    Compute::getChangeInPortfolioValueForAssetValueToBecomeIdealValue(weight, assetValues[i], portfolioValue)
                );
            },
            0.0
        );
}
}




JSArray<double> Solution::getResult(void) const noexcept
{
    if (this->assetCount <= 1)
        return JSArray<double>(this->assetCount, this->targetBuySell);

    if (Compute::isCurrencyEqual(this->targetBuySell, 0.0))
        return JSArray<double>(this->assetCount, 0.0);

    const bool isSellWholePortfolio = -this->targetBuySell > this->totalPortfolioValue || Compute::isCurrencyEqual(-this->targetBuySell, this->totalPortfolioValue);
    if (isSellWholePortfolio)
        return this->assetValues.map([](double val){return -val;});

    if (this->numberOfZeroWeights > 0)
        return this->zeroWeightHandler();

    if (this->canDoFullSingleOpRebalance())
        return this->targetWeights.map(
            [&](double weight, std::size_t i)
            {
                const double finalPortfolioValue = this->totalPortfolioValue + this->targetBuySell;
                return Compute::rebalanceAssetValue(weight, this->assetValues[i], finalPortfolioValue);
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
    if (this->isSell)
    {
        const double minSell = Compute::getMinTargetSellForFullPortfolioOnlySellRebalance(this->targetWeights, this->assetValues, this->totalPortfolioValue);
        const double targetSell = this->targetBuySell; // is a negative value
        return -targetSell > -minSell || Compute::isCurrencyEqual(targetSell, minSell);
    }

    const double minBuy = Compute::getMinTargetBuyForFullPortfolioOnlyBuyRebalance(this->targetWeights, this->assetValues, this->totalPortfolioValue);
    const double targetBuy = this->targetBuySell;
    return targetBuy > minBuy || Compute::isCurrencyEqual(targetBuy, minBuy);
}




std::size_t Solution::getMostOverWeightIndex(/* targetWeights, */ const JSArray<double>& currAssetValues, const double currPortfolioValue) const noexcept
{
    std::size_t overWeightIndex = 0;
    double overWeightValue = Compute::rebalanceAssetValue(this->targetWeights[0], currAssetValues[0], currPortfolioValue);
    for (std::size_t i = 1; i < this->assetCount; i += 1)
    {
        const double iValue = Compute::rebalanceAssetValue(this->targetWeights[i], currAssetValues[i], currPortfolioValue);
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
    double underWeightValue = Compute::rebalanceAssetValue(this->targetWeights[0], currAssetValues[0], currPortfolioValue);
    for (std::size_t i = 1; i < this->assetCount; i += 1)
    {
        const double iValue = Compute::rebalanceAssetValue(this->targetWeights[i], currAssetValues[i], currPortfolioValue);
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
            return Compute::getAssetValueChangeForEqRebalance(weight0, value0, weight1, value1, this->totalPortfolioValue);
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
            if (!Compute::isWeightZero(weight))
            {
                result[k] = i;
                k += 1;
            }
        }

        return result;
    })();
    const JSArray<double> targetNonZeroWeights = indexNumbersNonZeroWeight
        .map([&](std::size_t valIndex){return this->targetWeights[valIndex];});
    const JSArray<double> nonZeroWeightAssets = indexNumbersNonZeroWeight
        .map([&](std::size_t valIndex){return this->assetValues[valIndex];});
    const double sumOfNonZeroWeightAssets = nonZeroWeightAssets
        .reduce([](double sum, double val){return sum + val;}, 0.0);

    if (!this->isSell)
    {
        // is buy and just ignore buying anything with a 0% target weight
        return alignResultsToOriginalIndex(
            Solution(
                targetNonZeroWeights,
                nonZeroWeightAssets,
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
    const JSArray<std::size_t> indexNumbersZeroWeight = ([&]()
    {
        JSArray<std::size_t> result(this->numberOfZeroWeights);
        for (std::size_t i = 0, k = 0; i < this->assetCount; i += 1)
        {
            const double weight = this->targetWeights[i];
            if (Compute::isWeightZero(weight))
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
            const std::size_t originalChangePacket = this->targetBuySell / originalMaxIter;
            const std::size_t equivalentItersDoneBySellingAllZeroWeightAssets = sumOfZeroWeightAssetValues / originalChangePacket;
            return originalMaxIter - equivalentItersDoneBySellingAllZeroWeightAssets;
        })();
        const JSArray<double> changesToNonZeroAssetValues = Solution(
            targetNonZeroWeights,
            nonZeroWeightAssets,
            sumOfNonZeroWeightAssets,
            0,
            targetSell + sumOfZeroWeightAssetValues,
            newMaxIter
        ).getResult();
        JSArray<double> result = alignResultsToOriginalIndex(
            changesToNonZeroAssetValues,
            indexNumbersNonZeroWeight,
            this->assetCount
        );

        for (const std::size_t zeroIndex : indexNumbersZeroWeight)
        {
            result[zeroIndex] = -assetValues[zeroIndex];
        }

        return result;
    }

    if (canOnlySellSomeZeroWeightAssets)
    {

        JSArray<double> result = this->sellZeroWeightAssets(indexNumbersZeroWeight);
        for (std::size_t i = 0; i < this->numberOfZeroWeights; i += 1)
        {
            const double afterSaleValue = result[i];
            const double zeroWeightAssetValue = this->assetValues[indexNumbersZeroWeight[i]];
            // change result to represent how much to sell exactly, not the ending value after the sale
            result[i] = afterSaleValue - zeroWeightAssetValue;
        }

        return alignResultsToOriginalIndex(
            result,
            indexNumbersZeroWeight,
            targetWeights.size()
        );
    }

    // can exactly only sell the zero weight assets
    return ([&]() -> JSArray<double>
    {
        JSArray<double> result(this->assetCount, 0);
        for (const std::size_t zeroIndex : indexNumbersZeroWeight)
        {
            result[zeroIndex] = -assetValues[zeroIndex];
        }

        return result;
    })();
}





JSArray<double> Solution::sellZeroWeightAssets(const JSArray<std::size_t>& indexNumbersZeroWeight) const noexcept
{
    const double targetSell = this->targetBuySell; // negative value
    if (this->numberOfZeroWeights == 1) return {this->assetValues[indexNumbersZeroWeight[0]] + targetSell};

    JSArray<double> afterSalesZeroWeightAssetValues = indexNumbersZeroWeight
        .map([&](std::size_t zeroIndex){return this->assetValues[zeroIndex];});
    // used to be able to access in sorted fashion but still keep original order of values in result
    JSArray<double*> valuesSortedAccessorDescending = ([&]() -> JSArray<double*>
    {
        JSArray<double*> result(this->numberOfZeroWeights);
        for (std::size_t i = 0; i < this->numberOfZeroWeights; i += 1)
        {
            result[i] = &afterSalesZeroWeightAssetValues[i];
        }

        return result.sort([](const double* a, const double* b){return *a > *b;});
    })();
    double currLargestValue = *(valuesSortedAccessorDescending[0]);
    std::size_t groupSizeWithSameLargestValue = ([&]() -> std::size_t
    {
        std::size_t groupSize = 1;
        while (
            groupSize < this->numberOfZeroWeights &&
            Compute::isCurrencyEqual(
                currLargestValue,
                *(valuesSortedAccessorDescending[groupSize])
            )
        ){groupSize += 1;}

        return groupSize;
    })();
    double remainingSell = targetSell; // remember, this value is negative!!!!

    while (groupSizeWithSameLargestValue < this->numberOfZeroWeights && !Compute::isCurrencyEqual(remainingSell, 0.0))
    {
        const std::size_t startIndexOfNextGroup = groupSizeWithSameLargestValue;
        const double nextGroupUniformValue = *(valuesSortedAccessorDescending[startIndexOfNextGroup]);
        const double idealSell = -(currLargestValue - nextGroupUniformValue) * groupSizeWithSameLargestValue;
        const double actualSell = (-remainingSell > -idealSell) ? idealSell : remainingSell;

        remainingSell -= actualSell;
        currLargestValue += actualSell / groupSizeWithSameLargestValue;
        // maybe the size of the curr next group is greater than 1, this ensures we also get those values to.
        while (
            groupSizeWithSameLargestValue < this->numberOfZeroWeights &&
            Compute::isCurrencyEqual(
                *(valuesSortedAccessorDescending[groupSizeWithSameLargestValue]),
                currLargestValue
            )
        ){groupSizeWithSameLargestValue += 1;}
    }
    currLargestValue += remainingSell / groupSizeWithSameLargestValue;

    for (std::size_t i = 0; i < groupSizeWithSameLargestValue; i += 1)
    {
        *(valuesSortedAccessorDescending[i]) = currLargestValue;
    }

    return afterSalesZeroWeightAssetValues;
}
