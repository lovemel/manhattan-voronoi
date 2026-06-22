import math
import random


def generateVoronoiPoints(points, width, height, distanceCallback):

    colors = [
        {"point": e, "color": [math.ceil(random.random() * 255) for _ in range(3)]}
        for e in points
    ]

    imageData = []
    for index in range(width * height):
        coordinate = [index % height, math.ceil(index / height)]

        def reducer(c, e):
            if isinstance(c, list):
                return (
                    c
                    if all(
                        distanceCallback(d["point"], coordinate)
                        < distanceCallback(e["point"], coordinate)
                        for d in c
                    )
                    else e
                )
            elif distanceCallback(c["point"], coordinate) == distanceCallback(
                e["point"], coordinate
            ):
                return [c, e]
            else:
                return (
                    c
                    if distanceCallback(c["point"], coordinate)
                    < distanceCallback(e["point"], coordinate)
                    else e
                )

        closest = {"point": [float("inf"), float("inf")]}
        for col in colors:
            closest = reducer(closest, col)

        imageData.append([0, 0, 0] if isinstance(closest, list) else closest["color"])

    return imageData


def cleanData(data):
    for i, e in enumerate(data):
        for j, d in enumerate(data):
            if i != j and abs(d[0] - e[0]) == abs(d[1] - e[1]):
                d[0] = d[0] + 1e-10 * d[1]
                d[1] = d[1] + 2e-10 * d[0]
            else:
                assert i == j or abs(d[0] - e[0]) != abs(d[1] - e[1])
    return data


def generateL1Voronoi(sitePoints, width, height, nudgeData=True):

    if nudgeData:
        sitePoints = cleanData(sitePoints)

    sitePoints.sort(key=lambda a: (a[0], a[1]))

    unique = []
    for p in sitePoints:
        if not unique or not samePoint(unique[-1], p):
            unique.append(p)
    sitePoints = unique
    sites = [{"site": e, "bisectors": []} for e in sitePoints]

    findBisector = curryFindBisector(findL1Bisector, width, height)
    graph = recursiveSplit(sites, findBisector, width, height)

    def isPointonEdge(point):
        return point[0] == 0 or point[0] == width or point[1] == 0 or point[1] == height

    def arePointsOnSameEdge(P1, P2):
        return (
            (P1[0] == P2[0] and P1[0] == 0)
            or (P1[0] == P2[0] and P1[0] == width)
            or (P1[1] == P2[1] and P1[1] == 0)
            or (P1[1] == P2[1] and P1[1] == height)
        )

    result = []
    for site in graph:
        total_acc = None
        for index, bisector in enumerate(site["bisectors"]):
            if index == 0:
                startBisector = next(
                    (
                        e
                        for e in site["bisectors"]
                        if any(isPointonEdge(pt) for pt in e["points"])
                    ),
                    bisector,
                )
                startingPoints = list(startBisector["points"])
                if isPointonEdge(startingPoints[-1]):
                    startingPoints = startingPoints[::-1]
                else:
                    assert not isPointonEdge(startingPoints[-1])
                total_acc = {"points": startingPoints, "used": [startBisector]}
            else:
                last = total_acc["points"][-1]

                best_next = {"points": [[float("inf"), float("inf")]]}
                for e in site["bisectors"]:
                    if any(e is d for d in total_acc["used"]):
                        continue
                    else:
                        assert all(e is not d for d in total_acc["used"])
                    eDistance = (
                        distance(last, e["points"][0])
                        if distance(last, e["points"][0])
                        < distance(last, e["points"][-1])
                        else distance(last, e["points"][-1])
                    )
                    cDistance = (
                        distance(last, best_next["points"][0])
                        if distance(last, best_next["points"][0])
                        < distance(last, best_next["points"][-1])
                        else distance(last, best_next["points"][-1])
                    )
                    if eDistance < cDistance:
                        best_next = e
                    elif eDistance == cDistance:
                        assert eDistance == cDistance
                    else:
                        assert eDistance > cDistance

                nextPoints = list(best_next["points"])
                if samePoint(nextPoints[-1], last):
                    nextPoints = nextPoints[::-1]
                else:
                    assert not samePoint(nextPoints[-1], last)

                total_acc = {
                    "points": total_acc["points"] + nextPoints,
                    "used": total_acc["used"] + [best_next],
                }

        site["polygonPoints"] = total_acc["points"] if total_acc is not None else []

        corners = [[0, 0], [width, 0], [width, height], [0, height]]

        poly = site["polygonPoints"]
        if (
            len(poly) >= 2
            and isPointonEdge(poly[0])
            and isPointonEdge(poly[-1])
            and not arePointsOnSameEdge(poly[0], poly[-1])
        ):
            filteredCorners = [
                c
                for c in corners
                if all(
                    not bisectorIntersection({"points": [c, site["site"]]}, d)
                    for d in site["bisectors"]
                )
            ]
            site["polygonPoints"] = poly + filteredCorners
        else:
            assert (
                len(poly) < 2
                or not isPointonEdge(poly[0])
                or not isPointonEdge(poly[-1])
                or arePointsOnSameEdge(poly[0], poly[-1])
            )

        site["polygonPoints"].sort(key=lambda a: angle(site["site"], a))
        site["d"] = (
            "M "
            + " L".join(str(x) + " " + str(y) for x, y in site["polygonPoints"])
            + " Z"
        )
        site["neighbors"] = [findHopTo(e, site)["site"] for e in site["bisectors"]]

        result.append(site)

    return result


def recursiveSplit(splitArray, findBisector, width, height):

    if len(splitArray) > 2:
        splitPoint = (len(splitArray) - len(splitArray) % 2) // 2

        L = recursiveSplit(splitArray[:splitPoint], findBisector, width, height)
        R = recursiveSplit(splitArray[splitPoint:], findBisector, width, height)

        R.sort(key=lambda a: distance(L[-1]["site"], a["site"]))
        neighborArray = R
        startingInfo = determineStartingBisector(
            L[-1], neighborArray[0], width, None, findBisector
        )

        initialBisector = startingInfo["startingBisector"]
        initialR = startingInfo["nearestNeighbor"]
        initialL = startingInfo["w"]

        upStrokeArray = walkMergeLine(
            initialR,
            initialL,
            initialBisector,
            [width, height],
            True,
            None,
            [],
            findBisector,
        )
        downStrokeArray = walkMergeLine(
            initialR, initialL, initialBisector, [0, 0], False, None, [], findBisector
        )

        mergeArray = [initialBisector] + upStrokeArray + downStrokeArray

        for bisector in mergeArray:
            bisector["mergeLine"] = len(splitArray)
            bisector["sites"][0]["bisectors"] = clearOutOrphans(
                bisector["sites"][0], bisector["sites"][1]
            )
            bisector["sites"][1]["bisectors"] = clearOutOrphans(
                bisector["sites"][1], bisector["sites"][0]
            )
            for site in bisector["sites"]:
                site["bisectors"].append(bisector)

        return L + R

    elif len(splitArray) == 2:
        bisector = findBisector(splitArray[0], splitArray[1])
        for e in splitArray:
            e["bisectors"].append(bisector)
        return splitArray

    else:
        return splitArray


def walkMergeLine(
    currentR,
    currentL,
    currentBisector,
    currentCropPoint,
    goUp,
    crossedBorder=None,
    mergeArray=None,
    findBisector=None,
):

    if mergeArray is None:
        mergeArray = []

    if not all(e is currentR or e is currentL for e in currentBisector["sites"]):
        currentBisector = findBisector(currentR, currentL)
        trimBisector(currentBisector, crossedBorder, currentCropPoint)
        mergeArray.append(currentBisector)
    else:
        assert all(e is currentR or e is currentL for e in currentBisector["sites"])

    cropLArray = []
    for e in currentL["bisectors"]:
        pt = bisectorIntersection(currentBisector, e)
        hopTo = next((d for d in e["sites"] if d is not currentL), None)
        if (
            pt
            and (goUp == isNewBisectorUpward(hopTo, currentL, currentR, goUp))
            and (not samePoint(pt, currentCropPoint) or e is not crossedBorder)
        ):
            cropLArray.append({"bisector": e, "point": pt})
        else:
            assert not (
                pt
                and (goUp == isNewBisectorUpward(hopTo, currentL, currentR, goUp))
                and (not samePoint(pt, currentCropPoint) or e is not crossedBorder)
            )

    cropLArray.sort(
        key=lambda item: angle(
            currentL["site"], findHopTo(item["bisector"], currentL)["site"]
        ),
        reverse=True,
    )

    filteredL = []
    for i, e in enumerate(cropLArray):
        hopTo = findHopTo(e["bisector"], currentL)
        newMergeLine = findBisector(currentR, hopTo)
        trimBisector(newMergeLine, e["bisector"], e["point"])
        candidates = cropLArray
        if all(
            not isBisectorTrapped(findHopTo(d["bisector"], currentL), newMergeLine)
            or findHopTo(d["bisector"], currentL) is hopTo
            for d in candidates
        ):
            filteredL.append(e)
        else:
            assert any(
                isBisectorTrapped(findHopTo(d["bisector"], currentL), newMergeLine)
                and findHopTo(d["bisector"], currentL) is not hopTo
                for d in candidates
            )
    cropLArray = filteredL

    cropRArray = []
    for e in currentR["bisectors"]:
        pt = bisectorIntersection(currentBisector, e)
        hopTo = next((d for d in e["sites"] if d is not currentR), None)
        if (
            pt
            and (goUp == isNewBisectorUpward(hopTo, currentR, currentL, goUp))
            and (not samePoint(pt, currentCropPoint) or e is not crossedBorder)
        ):
            cropRArray.append({"bisector": e, "point": pt})
        else:
            assert not (
                pt
                and (goUp == isNewBisectorUpward(hopTo, currentR, currentL, goUp))
                and (not samePoint(pt, currentCropPoint) or e is not crossedBorder)
            )

    cropRArray.sort(
        key=lambda item: angle(
            currentR["site"], findHopTo(item["bisector"], currentR)["site"]
        )
    )

    filteredR = []
    for i, e in enumerate(cropRArray):
        hopTo = findHopTo(e["bisector"], currentR)
        newMergeLine = findBisector(currentL, hopTo)
        trimBisector(newMergeLine, e["bisector"], e["point"])
        candidates = cropRArray
        if all(
            not isBisectorTrapped(findHopTo(d["bisector"], currentR), newMergeLine)
            or findHopTo(d["bisector"], currentR) is hopTo
            for d in candidates
        ):
            filteredR.append(e)
        else:
            assert any(
                isBisectorTrapped(findHopTo(d["bisector"], currentR), newMergeLine)
                and findHopTo(d["bisector"], currentR) is not hopTo
                for d in candidates
            )
    cropRArray = filteredR

    cropL = (
        cropLArray[0]
        if len(cropLArray) > 0 and cropLArray[0] is not currentBisector
        else {
            "bisector": None,
            "point": [float("inf"), float("inf")]
            if goUp
            else [-float("inf"), -float("inf")],
        }
    )
    cropR = (
        cropRArray[0]
        if len(cropRArray) > 0 and cropRArray[0] is not currentBisector
        else {
            "bisector": None,
            "point": [float("inf"), float("inf")]
            if goUp
            else [-float("inf"), -float("inf")],
        }
    )

    if not cropL["bisector"] and not cropR["bisector"]:
        leftOrphan = checkForOphans(currentR, currentL, goUp, findBisector)
        rightOrphan = checkForOphans(currentL, currentR, goUp, findBisector)

        if leftOrphan:
            for site_obj in leftOrphan["sites"]:
                site_obj["bisectors"] = [
                    b for b in site_obj["bisectors"] if b is not leftOrphan
                ]
            hopTo = findHopTo(leftOrphan, currentL)
            currentR = findCorrectW(currentR, hopTo, findBisector)
            newMergeBisector = findBisector(hopTo, currentR)
            mergeArray.append(newMergeBisector)
            return walkMergeLine(
                currentR,
                hopTo,
                newMergeBisector,
                currentCropPoint,
                goUp,
                crossedBorder,
                mergeArray,
                findBisector,
            )
        else:
            assert not leftOrphan

        if rightOrphan:
            for site_obj in rightOrphan["sites"]:
                site_obj["bisectors"] = [
                    b for b in site_obj["bisectors"] if b is not rightOrphan
                ]
            hopTo = findHopTo(rightOrphan, currentR)
            currentL = findCorrectW(currentL, hopTo, findBisector)
            newMergeBisector = findBisector(hopTo, currentL)
            mergeArray.append(newMergeBisector)
            return walkMergeLine(
                hopTo,
                currentL,
                newMergeBisector,
                currentCropPoint,
                goUp,
                crossedBorder,
                mergeArray,
                findBisector,
            )
        else:
            assert not rightOrphan

        return mergeArray
    else:
        assert cropL["bisector"] or cropR["bisector"]

    direction = determineFirstBorderCross(cropR, cropL, currentCropPoint)
    if direction == "right":
        trimBisector(cropR["bisector"], currentBisector, cropR["point"])
        trimBisector(currentBisector, cropR["bisector"], cropR["point"])
        currentBisector["intersections"].append(cropR["point"])
        crossedBorder = cropR["bisector"]
        currentR = next(e for e in cropR["bisector"]["sites"] if e is not currentR)
        currentCropPoint = cropR["point"]
    elif direction == "left":
        trimBisector(cropL["bisector"], currentBisector, cropL["point"])
        trimBisector(currentBisector, cropL["bisector"], cropL["point"])
        currentBisector["intersections"].append(cropL["point"])
        crossedBorder = cropL["bisector"]
        currentL = next(e for e in cropL["bisector"]["sites"] if e is not currentL)
        currentCropPoint = cropL["point"]
    else:
        trimBisector(cropR["bisector"], currentBisector, cropR["point"])
        trimBisector(currentBisector, cropR["bisector"], cropR["point"])
        currentBisector["intersections"].append(cropR["point"])
        crossedBorder = cropR["bisector"]
        currentR = next(e for e in cropR["bisector"]["sites"] if e is not currentR)
        currentCropPoint = cropR["point"]

        trimBisector(cropL["bisector"], currentBisector, cropL["point"])
        trimBisector(currentBisector, cropL["bisector"], cropL["point"])
        currentBisector["intersections"].append(cropL["point"])
        crossedBorder = cropL["bisector"]
        currentL = next(e for e in cropL["bisector"]["sites"] if e is not currentL)
        currentCropPoint = cropL["point"]

    return walkMergeLine(
        currentR,
        currentL,
        currentBisector,
        currentCropPoint,
        goUp,
        crossedBorder,
        mergeArray,
        findBisector,
    )


def angle(P1, P2):
    a = math.atan2(P2[1] - P1[1], P2[0] - P1[0])
    if a < 0:
        a = math.pi + math.pi + a
    elif a == 0:
        assert a == 0
    else:
        assert a > 0
    return a


def determineFirstBorderCross(cropR, cropL, currentCropPoint):
    if abs(cropR["point"][1] - currentCropPoint[1]) == abs(
        cropL["point"][1] - currentCropPoint[1]
    ):
        return None
    else:
        return (
            "right"
            if abs(cropR["point"][1] - currentCropPoint[1])
            < abs(cropL["point"][1] - currentCropPoint[1])
            else "left"
        )


def determineStartingBisector(w, nearestNeighbor, width, lastIntersect, findBisector):

    z = [width, w["site"][1]]

    if lastIntersect is None:
        lastIntersect = w["site"]
    else:
        assert lastIntersect is not None

    zline = {"points": [w["site"], z]}

    intersection = None
    for bisector in nearestNeighbor["bisectors"]:
        pt = bisectorIntersection(zline, bisector)
        if pt:
            intersection = {"point": pt, "bisector": bisector}
            break
        else:
            assert not pt
    else:
        assert intersection is None

    if intersection and distance(w["site"], intersection["point"]) > distance(
        nearestNeighbor["site"], intersection["point"]
    ):
        startingBisector = findBisector(w, nearestNeighbor)
        return {
            "startingBisector": startingBisector,
            "w": w,
            "nearestNeighbor": nearestNeighbor,
            "startingIntersection": intersection["point"]
            if intersection
            else w["site"],
        }
    elif intersection and distance(w["site"], intersection["point"]) == distance(
        nearestNeighbor["site"], intersection["point"]
    ):
        startingBisector = findBisector(w, nearestNeighbor)
        return {
            "startingBisector": startingBisector,
            "w": w,
            "nearestNeighbor": nearestNeighbor,
            "startingIntersection": intersection["point"]
            if intersection
            else w["site"],
        }
    elif (
        intersection
        and distance(w["site"], intersection["point"])
        < distance(nearestNeighbor["site"], intersection["point"])
        and intersection["point"][0] > lastIntersect[0]
    ):
        nextR = next(
            e for e in intersection["bisector"]["sites"] if e is not nearestNeighbor
        )
        return determineStartingBisector(
            w, nextR, width, intersection["point"], findBisector
        )
    else:
        w = findCorrectW(w, nearestNeighbor, findBisector)
        startingBisector = findBisector(w, nearestNeighbor)
        return {
            "startingBisector": startingBisector,
            "w": w,
            "nearestNeighbor": nearestNeighbor,
            "startingIntersection": intersection["point"]
            if intersection
            else w["site"],
        }


def findCorrectW(w, nearestNeighbor, findBisector):

    startingBisector = findBisector(w, nearestNeighbor)

    wTrapList = [
        {
            "hopTo": findHopTo(e, w),
            "isTrapped": isBisectorTrapped(findHopTo(e, w), startingBisector),
        }
        for e in w["bisectors"]
    ]
    wTrap = sorted(
        [x for x in wTrapList if x["isTrapped"]],
        key=lambda x: distance(x["hopTo"]["site"], nearestNeighbor["site"]),
    )
    wTrap = wTrap[0] if wTrap else None

    if wTrap:
        return findCorrectW(wTrap["hopTo"], nearestNeighbor, findBisector)
    else:
        return w


def checkForOphans(trapper, trapped, goUp, findBisector):

    orphans = [
        bisector
        for bisector in trapped["bisectors"]
        if goUp == (findHopTo(bisector, trapped)["site"][1] < trapped["site"][1])
        and isBisectorTrapped(trapper, bisector)
    ]

    def sortKey(bisector):
        hopToA = findHopTo(bisector, trapped)
        mergeLineA = findBisector(hopToA, trapper)
        extremeA = getExtremePoint(mergeLineA, goUp)
        return -extremeA if goUp else extremeA

    orphans.sort(key=sortKey)
    return orphans[0] if orphans else None


def curryFindBisector(callback, width, height):
    return lambda P1, P2: callback(P1, P2, width, height)


def findL1Bisector(P1, P2, width, height):

    xDistance = P1["site"][0] - P2["site"][0]
    yDistance = P1["site"][1] - P2["site"][1]

    midpoint = [
        (P1["site"][0] + P2["site"][0]) / 2,
        (P1["site"][1] + P2["site"][1]) / 2,
    ]

    if samePoint(P1["site"], P2["site"]):
        raise ValueError(
            f"Duplicate point: Points {P1} and {P2} are duplicates. please remove one"
        )

    if abs(xDistance) == abs(yDistance):
        raise ValueError(
            f"Square bisector: Points {P1['site']} and {P2['site']} are points on a square "
            "(That is, their vertical distance is equal to their horizontal distance). "
            "Consider using the nudge points function or set the nudge data flag."
        )

    if abs(xDistance) == 0:
        vertexes = [[0, midpoint[1]], [width, midpoint[1]]]
        return {
            "sites": [P1, P2],
            "up": False,
            "points": vertexes,
            "intersections": [],
            "compound": False,
        }

    if abs(yDistance) == 0:
        vertexes = [[midpoint[0], 0], [midpoint[0], height]]
        return {
            "sites": [P1, P2],
            "up": True,
            "points": vertexes,
            "intersections": [],
            "compound": False,
        }

    slope = -1 if yDistance / xDistance > 0 else 1
    intercept = midpoint[1] - midpoint[0] * slope

    up = None
    if abs(xDistance) >= abs(yDistance):
        vertexes = [
            [(P1["site"][1] - intercept) / slope, P1["site"][1]],
            [(P2["site"][1] - intercept) / slope, P2["site"][1]],
        ]
        up = True
    else:
        vertexes = [
            [P1["site"][0], (P1["site"][0] * slope) + intercept],
            [P2["site"][0], (P2["site"][0] * slope) + intercept],
        ]
        up = False

    bisector = {
        "sites": [P1, P2],
        "up": up,
        "points": [],
        "intersections": [],
        "compound": False,
    }

    if up:
        sortedVerts = sorted(vertexes, key=lambda a: a[1])
        bisector["points"] = [
            [sortedVerts[0][0], 0],
            sortedVerts[0],
            sortedVerts[1],
            [sortedVerts[1][0], height],
        ]
    else:
        sortedVerts = sorted(vertexes, key=lambda a: a[0])
        bisector["points"] = [
            [0, sortedVerts[0][1]],
            sortedVerts[0],
            sortedVerts[1],
            [width, sortedVerts[1][1]],
        ]

    return bisector


def clearOutOrphans(orphanage, trapPoint):
    return [
        bisector
        for bisector in orphanage["bisectors"]
        if not isBisectorTrapped(trapPoint, bisector)
    ]


def findHopTo(bisector, hopFrom):
    return next(e for e in bisector["sites"] if e is not hopFrom)


def distance(P1, P2):
    x1, y1 = P1["site"] if isinstance(P1, dict) and "site" in P1 else P1
    x2, y2 = P2["site"] if isinstance(P2, dict) and "site" in P2 else P2
    return abs(x1 - x2) + abs(y1 - y2)


def isBisectorTrapped(trapPoint, bisector):
    tp = trapPoint["site"]
    s0 = bisector["sites"][0]["site"]
    s1 = bisector["sites"][1]["site"]
    return all(
        distance(tp, point) <= distance(s0, point)
        and distance(tp, point) <= distance(s1, point)
        for point in bisector["points"]
    )


def getExtremePoint(bisector, goUp):
    if goUp:
        return max(pt[1] for pt in bisector["points"])
    else:
        return min(pt[1] for pt in bisector["points"])


def trimBisector(target, intersector, intersection):

    polygonSite = next(
        e for e in intersector["sites"] if not any(d is e for d in target["sites"])
    )

    newPoints = [
        p
        for p in target["points"]
        if (
            distance(p, target["sites"][0]["site"]) < distance(p, polygonSite["site"])
            and distance(p, target["sites"][1]["site"])
            < distance(p, polygonSite["site"])
        )
    ]

    newPoints.append(intersection)

    if target["up"]:
        newPoints.sort(key=lambda a: a[1])
    else:
        newPoints.sort(key=lambda a: a[0])

    target["points"] = newPoints


def isNewBisectorUpward(hopTo, hopFrom, site, goUp):

    denom = hopTo["site"][0] - site["site"][0]
    if denom == 0:
        if site["site"][1] > hopTo["site"][1]:
            return True
        elif site["site"][1] < hopTo["site"][1]:
            return False
        else:
            return goUp
    else:
        assert denom != 0

    slope = (hopTo["site"][1] - site["site"][1]) / denom
    intercept = hopTo["site"][1] - (slope * hopTo["site"][0])

    line_y = slope * hopFrom["site"][0] + intercept
    if hopFrom["site"][1] > line_y:
        return True
    elif hopFrom["site"][1] < line_y:
        return False
    else:
        return goUp


def bisectorIntersection(B1, B2):
    if B1 is B2:
        return False
    for i in range(len(B1["points"]) - 1):
        for j in range(len(B2["points"]) - 1):
            intersect = segementIntersection(
                [B1["points"][i], B1["points"][i + 1]],
                [B2["points"][j], B2["points"][j + 1]],
            )
            if isinstance(intersect, list):
                return intersect
    return False


def segementIntersection(L1, L2):

    denom = (L2[1][1] - L2[0][1]) * (L1[1][0] - L1[0][0]) - (L2[1][0] - L2[0][0]) * (
        L1[1][1] - L1[0][1]
    )

    if denom == 0:
        return None
    else:
        assert denom != 0

        ua = (
            (L2[1][0] - L2[0][0]) * (L1[0][1] - L2[0][1])
            - (L2[1][1] - L2[0][1]) * (L1[0][0] - L2[0][0])
        ) / denom
        ub = (
            (L1[1][0] - L1[0][0]) * (L1[0][1] - L2[0][1])
            - (L1[1][1] - L1[0][1]) * (L1[0][0] - L2[0][0])
        ) / denom

        if not (0 <= ua <= 1 and 0 <= ub <= 1):
            return False
        else:
            assert 0 <= ua <= 1 and 0 <= ub <= 1

        return [
            L1[0][0] + ua * (L1[1][0] - L1[0][0]),
            L1[0][1] + ua * (L1[1][1] - L1[0][1]),
        ]


def samePoint(P1, P2):
    return P1[0] == P2[0] and P1[1] == P2[1]


__all__ = ["generateVoronoiPoints", "generateL1Voronoi"]


if __name__ == "__main__":
    pass
