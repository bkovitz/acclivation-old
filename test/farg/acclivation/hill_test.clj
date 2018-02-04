(ns farg.acclivation.hill-test
  (:refer-clojure :exclude [rand rand-int cond memoize])
  (:require [better-cond.core :refer [cond]]
            [clojure.tools.trace :refer [deftrace] :as trace]
            [clojure.pprint :refer [pprint]]
            [clojure.test :refer :all]
            [farg.acclivation.hill :refer :all]
            [farg.util :as util :refer [dd dde]]
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
  "Simple hill-shaped fitness function. Peak at origin."
  [xx]
  (->> xx
       (map #(- (Math/abs %)))
       (reduce +)))

(deftest test-hill-step
  (is (#{[-0.5 -1.0] [-1.0 -0.5]} (hill-step simple-f 0.5 [-1.0 -1.0])))
  (is (#{[-1.0 0.0] [-0.5 -0.5]} (hill-step simple-f 0.5 [-1.0 -0.5])))
  (is (= (hill-step simple-f 0.5 [0.0 0.0])
         [0.0 0.0])))

(deftest test-hill-climb
  (is (= (hill-climb simple-f 0.2 [-1.0 -1.0])
         ;fitness  n-steps  start xx    final xx
         [0.0      10       [-1.0 -1.0] [0.0 0.0]])))

(deftest test-deterministic-random-xxs
  (let [expect [[0.8 1.0] [-0.6 -0.8] [-0.9 -0.4] [-1.0 0.4]]]
    (is (= (take 4 (deterministic-random-xxs 0.1 2))
           expect))
    (is (= (take 4 (deterministic-random-xxs 0.1 2))
           expect))))

(deftest test-run-climbers
  (let [expect [[0.0 21 [0.8 1.0 -0.2 -0.1] [0.0 0.0 0.0 0.0]]
                [0.0 21 [-0.6 -0.8 -0.6 -0.1] [0.0 0.0 0.0 0.0]]
                [0.0 21 [-0.9 -0.4 -0.1 -0.7] [0.0 0.0 0.0 0.0]]
                [0.0 34 [-1.0 0.4 -1.0 -1.0] [0.0 0.0 0.0 0.0]]
                [0.0 21 [1.0 0.1 0.8 0.2] [0.0 0.0 0.0 0.0]]]]
    (is (= expect (run-climbers simple-f 0.1 4 5)))
    (is (= expect (run-climbers simple-f 0.1 4 5)))))
