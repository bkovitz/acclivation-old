(ns farg.acclivation.sa-test
  (:refer-clojure :exclude [rand rand-int cond memoize])
  (:require [better-cond.core :refer [cond]]
            [clojure.edn :as edn]
            [clojure.tools.trace :refer [deftrace] :as trace]
            [clojure.pprint :refer [pprint]]
            [clojure.test :refer :all]
            [com.rpl.specter :refer :all]
            [farg.acclivation.sa :refer :all]
            [farg.util :as util :refer [dd dde]]
            [farg.with-state :refer [with-state]]
            [farg.x.navs :refer :all]
            [ubergraph.core :as uber]))

(defn graph->map [g]
  {:nodes (into {} (->> (uber/nodes g)
                        (map #(vector % (uber/attrs g %)))))
   :edges (apply hash-set (->> (uber/edges g)
                               (map #(vector (uber/src %)
                                             (uber/dest %)
                                             (uber/attr g % :weight)))))})

(defn mrounded [g]
  (transform [(walker float?)] util/mround (graph->map g)))

(deftest test-sa
  (let [g0 (uber/digraph [:a :b -1.0])
        g1 (spread-activation g0 {:a 0.5} :iterations 1 :decay 0.5)
        g2 (spread-activation g1
             {:a (uber/attr g1 :a :a), :b (uber/attr g1 :b :a)}
             :iterations 1 :decay 0.5)
        g2' (spread-activation g0 {:a 0.5} :iterations 2 :decay 0.5)]
    (is (= g2 g2'))
    (is (= {:nodes {:a {:a 0.5}, :b {:a -0.268}}, :edges #{[:a :b -1.0]}}
           (mrounded g1)))
    (is (= {:nodes {:a {:a 0.5}, :b {:a -0.515}}, :edges #{[:a :b -1.0]}}
           (mrounded g2)))))
