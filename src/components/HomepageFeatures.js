import React from 'react';
import clsx from 'clsx';
import styles from './HomepageFeatures.module.css';

const FeatureList = [
  {
    title: 'Focus on What Matters',
    Svg: require('../../static/img/undraw_docusaurus_mountain.svg').default,
    description: (
      <>
        fenics-optimize provides reusable and straightforward coding
        for topology optimization problem. 
      </>
    ),
  },
  {
    title: 'Powerful solvers',
    Svg: require('../../static/img/undraw_docusaurus_tree.svg').default,
    description: (
      <>
        fenics-optimize consists of the high-performance optimization solvers
        such as the IPOPT-HSL and MMA solver.
      </>
    ),
  },
];

function Feature({Svg, title, description}) {
  return (
    <div className={clsx('col')}>
      <div className="text--center padding-horiz--md">
        <h3>{title}</h3>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
