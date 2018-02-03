(ns farg.acclivation.hill-test
  (:refer-clojure :exclude [rand rand-int cond memoize])
  (:require [better-cond.core :refer [cond]]
            [clojure.test :refer :all]
            [farg.acclivation.hill :refer :all]
            [farg.with-state :refer [with-state]]))

(deftest test-lazy-deterministic-shuffled-indices
  (is (= (lazy-deterministic-shuffled-indices 1 8)
         [5 0 3 1 2 7 4 6]))
  (is (= (lazy-deterministic-shuffled-indices 1 8)
         (lazy-deterministic-shuffled-indices 1 8)))
  (is (not= (lazy-deterministic-shuffled-indices 1 8)
            (lazy-deterministic-shuffled-indices 2 8)))
  (is (= (take 8 (lazy-deterministic-shuffled-indices 1 10000000))
         [9548985 5764588 641847 4970313 6064254 7814904 4504434 4906606])))

(deftest test-make-flatten-coords
  (let [cf (make-flatten-coords [-1.0 0.0 1.0] 2)]
    (is (= (map cf (range 9))
        [[-1.0 -1.0]
         [0.0 -1.0]
         [1.0 -1.0]
         [-1.0 0.0]
         [0.0 0.0]
         [1.0 0.0]
         [-1.0 1.0]
         [0.0 1.0]
         [1.0 1.0]]))))

(deftest test-all-steps
  (is (= (all-steps 0.5 [0.0 0.0])
         [[-0.5 -0.0] [0.5 -0.0] [0.0 -0.5] [0.0 0.5]]))
  (is (= (all-steps 0.5 [0.0 -1.0])
         [[-0.5 -1.0] [0.5 -1.0] [0.0 -0.5]])))

(defn simple-f
  "Simple hill-shaped fitness function. Peak at [0 0]."
  [[x0 x1]] (+ (- (Math/abs x0)) (- (Math/abs x1))))

(deftest test-hill-step
  (is (= (hill-step simple-f 0.5 [-1.0 -1.0])
         [-1.0 -0.5]))
  (is (= (hill-step simple-f 0.5 [-1.0 -0.5])
         [-1.0 0.0]))
  (is (= (hill-step simple-f 0.5 [0.0 0.0])
         [0.0 0.0])))
