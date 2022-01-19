import React from 'react';
import clsx from 'clsx';
import styles from './HomepageFeatures.module.css';

const FeatureList = [
  {
    title: 'Topology optimization in Python',
    Svg: require('../../static/img/undraw_docusaurus_mountain.svg').default,
    description: (
      <>
        Python is the most popular, easy-to-use, high-performance program language in the sciense.
        fenics-optimize provides toolbox for the general-purpose topology optimization in the python.
      </>
    ),
  },
  {
    title: 'Powerful solvers',
    Svg: require('../../static/img/undraw_docusaurus_tree.svg').default,
    description: (
      <>
        fenics-optimize contains high-performance optimization solvers such as the IPOPT-HSL and MMA solver.
      </>
    ),
  },
  {
    title: 'Parallel computing',
    Svg: require('../../static/img/undraw_docusaurus_tree.svg').default,
    description: (
      <>
        fenics-optimize has been designed for parallel processing through MPI.
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
