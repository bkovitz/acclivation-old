(ns farg.acclivation.hill
  "Hill-climbing, to measure acclivity of a function.

  All coordinates are in the range [-1.0, +1.0].

  Naming convention: xx is a vector of coordinates."
  (:refer-clojure :exclude [rand rand-int cond memoize])
  (:require [better-cond.core :refer [cond]]
            [clojure.tools.trace :refer [deftrace] :as trace]
            [clojure.pprint :refer [pprint]]
            [clojure.math.combinatorics :as combo]
            [clojure.math.numeric-tower :as math]
            [clojure.java.io :refer [file writer]]
            [clj-time.local :as ltime]
            [farg.util :refer [dd dde choose choose-one with-rng-seed mround
                               defopts with-*out*]
             :as util]
            [farg.with-state :refer [with-state]]))

(defn lazy-deterministic-shuffled-indices [seed ub]
  (let [rng (java.util.Random. seed)]
    (letfn [(next-int [] (.nextInt rng ub))
            (lazy [so-far]
              (lazy-seq
                (when (< (count so-far) ub)
                  (let [i (->> (repeatedly next-int)
                               (filter #(not (contains? so-far %)))
                               (first))]
                    (cons i (lazy (conj so-far i)))))))]
      (lazy #{}))))

(defn make-flatten-coords
  "Returns a function that maps an integer index to a multidimensional
  coordinate.

  row: A vector of all the numbers in a single dimension.

  dimension: The number of dimensions in the coordinates."
  [row dimension]
  (let [row-size (count row)]
    (fn [index]
      (loop [index index, coords [], n dimension]
        (cond
          (zero? n) coords
          :let [i (rem index row-size)
                index (quot index row-size)]
          (recur index (conj coords (get row i)) (dec n)))))))

(defn in-range? [xx]
  (every? #(<= -1.0 % 1.0) xx))

(defn all-steps
  "Returns all coordinates that vary from 'xx' by 'step' in one dimension."
  [step xx]
  (->> (for [n (range (count xx))]
         [(update xx n - step) (update xx n + step)])
       (apply concat)
       (filter in-range?)))

(defn hill-step
  "Returns point with highest f within one step of xx. If none have higher
  f than xx, returns xx."
  [f step xx]
  (let [new-xx (apply max-key f (all-steps step xx))]
    (if (> (f new-xx) (f xx)) new-xx xx)))
