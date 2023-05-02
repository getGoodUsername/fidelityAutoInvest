#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <numeric>

#include <iostream>


std::vector<double> rebalance(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue);
std::vector<double> changeInTotalPortfolioValueNeededForCurrentAssetValuesToBecomeTheTargetValues(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue);
bool canDoFullSingleOpRebalance(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue, double targetBuySell);
std::size_t getMostOverWeightIndex(const std::vector<double>& targetWeights, const std::vector<double>& currAssetValues, double currPortfolioValue);
std::size_t getMostUnderWeightIndex(const std::vector<double>& targetWeights, const std::vector<double>& currAssetValues, double currPortfolioValue);
std::vector<double> singleOperationTransaction(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double targetBuySell, std::size_t maxIter);
std::vector<double> fast_singleOperationTransaction(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double targetBuySell, std::size_t maxIter);
std::vector<double> zeroWeightHandler(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue, double targetBuySell, std::size_t maxIter);
std::vector<double> sellZeroWeightAsset(const std::vector<double>& zeroWeightAssetValues, const std::vector<std::size_t>& zeroWeightOriginalIndexNumbers, double targetSell);
bool floatingPointEquality(double a, double b, double epsilon);
bool currencyEquality(double a, double b);
bool isWeightZero(double weight);

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

    std::cout << singleOperationTransaction(targetWeights, assetValues, targetBuySell, abs(targetBuySell)) << std::endl;

    return 0;
}





bool floatingPointEquality(double a, double b, double epsilon)
{
    return std::abs(a - b) <= epsilon;
}



bool currencyEquality(double a, double b)
{
    // min value in a dollar is 0.01 (aka a penny). epsilon of a tenth of a penny should be good enough.
    constexpr double epsilon = 0.001;
    return floatingPointEquality(a, b, epsilon);
}



bool isWeightZero(double weight)
{
    // I only really want weights with three decimal digits max. (ex. 5.535%). so if weight is less than or eq to 0.0001% I want to treat as zero
    constexpr double epsilon = 1e-6;
    return floatingPointEquality(weight, 0.0, epsilon);
}



/**
 * @brief weights can NOT equal zero.
 * 
 * @param targetWeights 
 * @param assetValues 
 * @param totalPortfolioValue 
 * @return std::vector<double> 
 */
std::vector<double> changeInTotalPortfolioValueNeededForCurrentAssetValuesToBecomeTheTargetValues(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue)
{
    std::vector<double> result(targetWeights.size());
    for (std::size_t i = 0; i < targetWeights.size(); i += 1)
    {
        const double weight = targetWeights[i];
        const double value = assetValues[i];
        result[i] = (value / weight) - totalPortfolioValue;
    }

    return result;
}





bool canDoFullSingleOpRebalance(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue, double targetBuySell)
{
    const bool isSell = targetBuySell < 0;
    const std::vector<double> change = changeInTotalPortfolioValueNeededForCurrentAssetValuesToBecomeTheTargetValues(targetWeights, assetValues, totalPortfolioValue);

    if (isSell)
    {
        const double minPortfolioSellForFullOnlySellRebalance = *std::min_element(change.begin(), change.end());
        return targetBuySell < minPortfolioSellForFullOnlySellRebalance || currencyEquality(targetBuySell, minPortfolioSellForFullOnlySellRebalance);
    }

    const double minPortfolioBuyForFullOnlyBuyRebalance = *std::max_element(change.begin(), change.end());
    return targetBuySell > minPortfolioBuyForFullOnlyBuyRebalance || currencyEquality(targetBuySell, minPortfolioBuyForFullOnlyBuyRebalance);
}





std::vector<double> rebalance(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue)
{
    std::vector<double> result(targetWeights.size());

    for (std::size_t i = 0; i < targetWeights.size(); i += 1)
    {
        const double idealValue = totalPortfolioValue * targetWeights[i];
        result[i] = idealValue - assetValues[i];
    }

    return result;
}





std::size_t getMostOverWeightIndex(const std::vector<double>& targetWeights, const std::vector<double>& currAssetValues, double currPortfolioValue)
{
    const std::vector<double> currRebalance = rebalance(targetWeights, currAssetValues, currPortfolioValue);

    std::size_t overWeightIndex = 0;
    for (std::size_t i = 1; i < currRebalance.size(); i += 1)
    {
        overWeightIndex = (currRebalance[i] < currRebalance[overWeightIndex]) ? i : overWeightIndex;
    }
    return overWeightIndex;
}





std::size_t getMostUnderWeightIndex(const std::vector<double>& targetWeights, const std::vector<double>& currAssetValues, double currPortfolioValue)
{
    const std::vector<double> currRebalance = rebalance(targetWeights, currAssetValues, currPortfolioValue);

    std::size_t underWeightIndex = 0;
    for (std::size_t i = 1; i < currRebalance.size(); i += 1)
    {
        underWeightIndex = (currRebalance[i] > currRebalance[underWeightIndex]) ? i : underWeightIndex;
    }
    return underWeightIndex;
}





std::vector<double> singleOperationTransaction(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double targetBuySell, std::size_t maxIter)
{
    const double totalPortfolioValue = ([&]() -> double
    {
        double sum = 0;
        for (std::size_t i = 0; i < assetValues.size(); i += 1)
        {
            sum += assetValues[i];
        }

        return sum;
    })();

    const bool thereAreZeroTargetWeights = ([&]()
    {
        bool result = false;
        for (std::size_t i = 0; i < targetWeights.size() && !result; i += 1)
        {
            result |= isWeightZero(targetWeights[i]);
        }

        return false;
    })();

    if (thereAreZeroTargetWeights)
        return zeroWeightHandler(targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter);

    if (canDoFullSingleOpRebalance(targetWeights, assetValues, totalPortfolioValue, targetBuySell))
        return rebalance(targetWeights, assetValues, totalPortfolioValue);


    const bool isSell = targetBuySell < 0;
    std::vector<double> endAssetValues = assetValues;
    const double changePacket = targetBuySell / maxIter;
    const auto assetSelector = (isSell) ? getMostOverWeightIndex : getMostUnderWeightIndex;
    const std::size_t initTargetAssetIndex = assetSelector(targetWeights, assetValues, totalPortfolioValue);
    const double initChange = ([&]() -> double
    {
        const double weight0 = targetWeights[initTargetAssetIndex];
        const double value0 = assetValues[initTargetAssetIndex];
        std::vector<double> assetValueChangesForEqRebalanceVal(targetWeights.size());
        const auto getValue = [&](std::size_t i) -> double
        {
            const double weight1 = targetWeights[i];
            const double value1 = assetValues[i];
            return (totalPortfolioValue * (weight1 - weight0) + value0 - value1) / (weight0 - weight1 - 1);
        };

        for (std::size_t i = 0; i < initTargetAssetIndex; i += 1)
        {
            assetValueChangesForEqRebalanceVal[i] = getValue(i);
        }

        // do this for convenience, actual value would be 0 but thats
        // an irrelevant result since want to find the value such that
        // it would make this most under/over weight index tied with
        // in place with another index.
        assetValueChangesForEqRebalanceVal[initTargetAssetIndex] = targetBuySell;

        for (std::size_t i = initTargetAssetIndex + 1; i < targetWeights.size(); i += 1)
        {
            assetValueChangesForEqRebalanceVal[i] = getValue(i);
        }

        return (isSell) ?
            *std::max_element(assetValueChangesForEqRebalanceVal.begin(), assetValueChangesForEqRebalanceVal.end()) :
            *std::min_element(assetValueChangesForEqRebalanceVal.begin(), assetValueChangesForEqRebalanceVal.end());
    })();
    double currPortfolioValue = totalPortfolioValue + initChange;
    const std::size_t remainingIter = (targetBuySell - initChange) / changePacket;

    endAssetValues[initTargetAssetIndex] += initChange;
    for (std::size_t i = 0; i < remainingIter; i += 1, currPortfolioValue += changePacket)
    {
        endAssetValues[assetSelector(targetWeights, endAssetValues, currPortfolioValue)] += changePacket;
    }

    const double targetFinalPortfolioValue = totalPortfolioValue + targetBuySell;
    endAssetValues[assetSelector(targetWeights, endAssetValues, currPortfolioValue)] += targetFinalPortfolioValue - currPortfolioValue;

    return ([&]() -> std::vector<double>
    {
        std::vector<double>& diff = endAssetValues;
        for (std::size_t i = 0; i < diff.size(); i += 1)
        {
            diff[i] -= assetValues[i];
        }

        return diff;
    })();
}




/**
 * @brief SUM(zeroWeightAssetValues) > -targetSell must be true.
 * 
 * @param zeroWeightAssetValues 
 * @param zeroWeightOriginalIndexNumbers 
 * @param targetSell 
 * @return std::vector<double> which is the new values (NOT THE DIFF) of the assets after selling.
 */
std::vector<double> sellZeroWeightAsset(const std::vector<double>& zeroWeightAssetValues, const std::vector<std::size_t>& zeroWeightOriginalIndexNumbers, double targetSell)
{
    if (zeroWeightAssetValues.size() == 1) return {zeroWeightAssetValues[0] + targetSell};

    // I need to be able to sort and still know what the element's
    // original index was in order to put them back into their
    // original order.
    struct ValueIndexPair
    {
        double value;
        double index;
    };

    std::vector<ValueIndexPair> zeroWeightAssetsSortedValuesDescending = ([&]() -> std::vector<ValueIndexPair>
    {
        std::vector<ValueIndexPair> result(zeroWeightAssetValues.size());
        for (std::size_t i = 0; i < zeroWeightAssetValues.size(); i += 1)
        {
            result[i].value = zeroWeightAssetValues[i];
            result[i].index = zeroWeightOriginalIndexNumbers[i];
        }

        std::sort(result.begin(), result.end(), [](const ValueIndexPair& a, const ValueIndexPair& b) -> bool {return a.value > b.value;});
        return result;
    })();
    std::size_t groupWithSameLargestValueSize = ([&]() -> std::size_t
    {
        std::size_t groupSize = 1;
        const double largestValue = zeroWeightAssetsSortedValuesDescending[0].value;
        for (;
            groupSize < zeroWeightAssetsSortedValuesDescending.size() &&
                currencyEquality(largestValue, zeroWeightAssetsSortedValuesDescending[groupSize].value);
            groupSize += 1
        ){}

        return groupSize;
    })();
    double remainingSell = targetSell; // remember, this value is negative!!!!
    double currLargestValue = zeroWeightAssetsSortedValuesDescending[0].value;

    while (groupWithSameLargestValueSize < zeroWeightAssetsSortedValuesDescending.size() && !currencyEquality(remainingSell, 0))
    {
        // todo, have to check the size of the next group...    
        const double nextGroupUniformValue = zeroWeightAssetsSortedValuesDescending[groupWithSameLargestValueSize].value;
        const double targetChange = -(currLargestValue - nextGroupUniformValue) * groupWithSameLargestValueSize;
        const double actualChange = (-remainingSell > -targetChange || currencyEquality(remainingSell, targetChange)) ? targetChange : remainingSell;

        remainingSell -= actualChange;
        currLargestValue += actualChange / groupWithSameLargestValueSize;
        // maybe the size of the next group is greater than 1, this ensures we also get those values to.
        for (;
            groupWithSameLargestValueSize < zeroWeightAssetsSortedValuesDescending.size() &&
                currencyEquality(zeroWeightAssetsSortedValuesDescending[groupWithSameLargestValueSize].value, currLargestValue);
            groupWithSameLargestValueSize += 1
        ){}
    }

    currLargestValue += remainingSell / groupWithSameLargestValueSize;

    for (std::size_t i = 0; i < groupWithSameLargestValueSize; i += 1)
    {
        zeroWeightAssetsSortedValuesDescending[i].value = currLargestValue;
    }

    // revert assets to their original order. index in ascending order
    std::sort(
        zeroWeightAssetsSortedValuesDescending.begin(),
        zeroWeightAssetsSortedValuesDescending.end(),
        [](const ValueIndexPair& a, const ValueIndexPair& b) -> bool {return a.index < b.index;}
    );

    return ([&]() -> std::vector<double>
    {
        std::vector<double> result(zeroWeightAssetValues.size());
        for (std::size_t i = 0 ; i < zeroWeightAssetValues.size(); i += 1)
        {
            result[i] = zeroWeightAssetsSortedValuesDescending[i].value;
        }

        return result;
    })();
}





std::vector<double> zeroWeightHandler(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double totalPortfolioValue, double targetBuySell, std::size_t maxIter)
{
    if (-targetBuySell > totalPortfolioValue || currencyEquality(-targetBuySell, totalPortfolioValue)) // will only work if a sell and therefore a negative value
    {
        std::vector<double> fullSell(assetValues.size());
        std::transform(assetValues.begin(), assetValues.end(), fullSell.begin(), [](double value){return -value;});
        return fullSell;
    }

    const auto alignResultsToOriginalIndex = [](
        const std::vector<double>& subsetValues,
        const std::vector<std::size_t>& originalSubsetElementIndex,
        std::size_t superSetSize
    ) -> std::vector<double>
    {
        std::vector<double> aligned(superSetSize);
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

    const std::vector<std::size_t> indexNumbersNonZeroWeight = ([&]() -> std::vector<std::size_t>
    {
        std::vector<std::size_t> result;
        std::remove_copy_if(targetWeights.begin(), targetWeights.end(), std::back_inserter(result), isWeightZero);
        return result;
    })();
    const std::vector<double> targetNonZeroWeights = ([&]() -> std::vector<double>
    {
        std::vector<double> result(indexNumbersNonZeroWeight.size());
        for (std::size_t i = 0; i < indexNumbersNonZeroWeight.size(); i += 1)
        {
            result[i] = targetWeights[indexNumbersNonZeroWeight[i]];
        }

        return result;
    })();
    const std::vector<double> assetValuesWithTargetNonZeroWeights = ([&]() -> std::vector<double>
    {
        std::vector<double> result(indexNumbersNonZeroWeight.size());
        for (std::size_t i = 0; i < indexNumbersNonZeroWeight.size(); i += 1)
        {
            result[i] = assetValues[indexNumbersNonZeroWeight[i]];
        }

        return result;
    })();;

    // though target buy sell should be greater than 0 not equal...
    // we do this to be consistent with the negation of isSell
    const bool isBuy = targetBuySell > 0.0 || currencyEquality(targetBuySell, 0.0);

    if (isBuy)
    {
        return alignResultsToOriginalIndex(
            singleOperationTransaction(targetNonZeroWeights, assetValuesWithTargetNonZeroWeights, targetBuySell, maxIter),
            indexNumbersNonZeroWeight,
            targetWeights.size()
        );
    }
    // else is sell

    const double targetSell = targetBuySell;
    const std::size_t numberOfZeroWeights = targetWeights.size() - indexNumbersNonZeroWeight.size();
    const std::vector<std::size_t> indexNumbersZeroWeight = ([&]() -> std::vector<std::size_t>
    {
        std::vector<std::size_t> result(numberOfZeroWeights);
        for (std::size_t i = 0, k = 0; i < targetWeights.size(); i += 1)
        {
            if (isWeightZero(targetWeights[i]))
            {
                result[k] = i;
                k += 1;
            }
        }

        return result;
    })();
    const double sumOfZeroWeightAssetValues = ([&]() -> double
    {
        double sum = 0;
        for (const std::size_t zeroIndex : indexNumbersZeroWeight)
        {
            sum += assetValues[zeroIndex];
        }

        return sum;
    })();
    const bool canSellMoreThanJustAllTheZeroWeightAssets = -targetSell > sumOfZeroWeightAssetValues;
    const bool canOnlySellSomeZeroWeightAssets = sumOfZeroWeightAssetValues > -targetSell;

    if (canSellMoreThanJustAllTheZeroWeightAssets)
    {
        const std::size_t newMaxIter = ([&]() -> std::size_t
        {
            const std::size_t originalMaxIter = maxIter;
            const std::size_t originalChangePacket = totalPortfolioValue / originalMaxIter;
            const std::size_t equivalentItersDoneBySellingAllZeroWeightAssets = sumOfZeroWeightAssetValues / originalChangePacket;
            return originalMaxIter - equivalentItersDoneBySellingAllZeroWeightAssets;
        })();
        const std::vector<double> changesToNonZeroAssetValues = singleOperationTransaction(
            targetNonZeroWeights,
            assetValuesWithTargetNonZeroWeights,
            targetSell + sumOfZeroWeightAssetValues,
            newMaxIter
        );
        return ([&]() -> std::vector<double>
        {
            std::vector<double> result = alignResultsToOriginalIndex(changesToNonZeroAssetValues, indexNumbersNonZeroWeight, targetWeights.size());
            for (const std::size_t zeroIndex : indexNumbersZeroWeight)
            {
                result[zeroIndex] = -assetValues[zeroIndex];
            }

            return result;
        })();
    }

    if (canOnlySellSomeZeroWeightAssets)
    {
        const std::vector<double> zeroWeightAssetValues = ([&]() -> std::vector<double>
        {
            std::vector<double> result(numberOfZeroWeights);
            for (std::size_t i = 0; i < indexNumbersZeroWeight.size(); i += 1)
            {
                result[i] = assetValues[indexNumbersZeroWeight[i]];
            }

            return result;
        })();
        return ([&]() -> std::vector<double>
        {
            std::vector<double> result = sellZeroWeightAsset(zeroWeightAssetValues, indexNumbersZeroWeight, targetSell);
            for (std::size_t i = 0; i < numberOfZeroWeights; i += 1)
            {
                const double afterSaleValue = result[i];
                result[i] = afterSaleValue - zeroWeightAssetValues[i];
            }

            return alignResultsToOriginalIndex(
                result,
                indexNumbersZeroWeight,
                targetWeights.size()
            );
        })();
    }

    // can sell exactly only fully the zero weight assets
    return ([&]() -> std::vector<double>
    {
        std::vector<double> result(0);
        for (const std::size_t zeroIndex : indexNumbersZeroWeight)
        {
            result[zeroIndex] = -assetValues[zeroIndex];
        }

        return result;
    })();
}





std::vector<double> fast_singleOperationTransaction(const std::vector<double>& targetWeights, const std::vector<double>& assetValues, double targetBuySell, std::size_t maxIter)
{
    const double totalPortfolioValue = std::accumulate(assetValues.begin(), assetValues.end(), 0.0);
    const std::size_t numberOfZeroWeights = std::count_if(targetWeights.begin(), targetWeights.end(), isWeightZero);

    if (numberOfZeroWeights > 0)
        return zeroWeightHandler(targetWeights, assetValues, totalPortfolioValue, targetBuySell, maxIter);

    if (canDoFullSingleOpRebalance(targetWeights, assetValues, totalPortfolioValue, targetBuySell))
        return rebalance(targetWeights, assetValues, totalPortfolioValue);


    const bool isSell = targetBuySell < 0;
    std::vector<double> endAssetValues = assetValues;
    const double changePacket = targetBuySell / maxIter;
    const auto assetSelector = (isSell) ? getMostOverWeightIndex : getMostUnderWeightIndex;
    const std::size_t initTargetAssetIndex = assetSelector(targetWeights, assetValues, totalPortfolioValue);
    const double initChange = ([&]() -> double
    {
        const double weight0 = targetWeights[initTargetAssetIndex];
        const double value0 = assetValues[initTargetAssetIndex];
        std::vector<double> assetValueChangesForEqRebalanceVal(targetWeights.size());
        const auto getValue = [&](std::size_t i) -> double
        {
            const double weight1 = targetWeights[i];
            const double value1 = assetValues[i];
            return (totalPortfolioValue * (weight1 - weight0) + value0 - value1) / (weight0 - weight1 - 1);
        };

        for (std::size_t i = 0; i < initTargetAssetIndex; i += 1)
        {
            assetValueChangesForEqRebalanceVal[i] = getValue(i);
        }

        // do this for convenience, actual value would be 0 but thats
        // an irrelevant result since want to find the value such that
        // it would make this most under/over weight index tied with
        // in place with another index.
        assetValueChangesForEqRebalanceVal[initTargetAssetIndex] = targetBuySell;

        for (std::size_t i = initTargetAssetIndex + 1; i < targetWeights.size(); i += 1)
        {
            assetValueChangesForEqRebalanceVal[i] = getValue(i);
        }

        return (isSell) ?
            *std::max_element(assetValueChangesForEqRebalanceVal.begin(), assetValueChangesForEqRebalanceVal.end()) :
            *std::min_element(assetValueChangesForEqRebalanceVal.begin(), assetValueChangesForEqRebalanceVal.end());
    })();
    double currPortfolioValue = totalPortfolioValue + initChange;
    const std::size_t remainingIter = (targetBuySell - initChange) / changePacket;

    endAssetValues[initTargetAssetIndex] += initChange;
    for (std::size_t i = 0; i < remainingIter; i += 1, currPortfolioValue += changePacket)
    {
        endAssetValues[assetSelector(targetWeights, endAssetValues, currPortfolioValue)] += changePacket;
    }

    const double targetFinalPortfolioValue = totalPortfolioValue + targetBuySell;
    endAssetValues[assetSelector(targetWeights, endAssetValues, currPortfolioValue)] += targetFinalPortfolioValue - currPortfolioValue;

    return ([&]() -> std::vector<double>
    {
        std::vector<double>& diff = endAssetValues;
        for (std::size_t i = 0; i < diff.size(); i += 1)
        {
            diff[i] -= assetValues[i];
        }

        return diff;
    })();
}
