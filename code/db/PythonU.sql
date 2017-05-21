-- phpMyAdmin SQL Dump
-- version 4.3.8
-- http://www.phpmyadmin.net
--
-- Host: localhost
-- Generation Time: Feb 23, 2017 at 01:44 AM
-- Server version: 5.5.51-38.2
-- PHP Version: 5.4.31

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Database: `PythonU`
--

-- --------------------------------------------------------

--
-- Table structure for table `Courses`
--

CREATE TABLE IF NOT EXISTS `Courses` (
  `CourseID` int(10) unsigned NOT NULL,
  `Course_Name` varchar(150) COLLATE utf8_unicode_ci NOT NULL
) ENGINE=MyISAM AUTO_INCREMENT=3 DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci COMMENT='Courses in Python U.';

--
-- Dumping data for table `Courses`
--

INSERT INTO `Courses` (`CourseID`, `Course_Name`) VALUES
(1, 'Python 101'),
(2, 'Mathematics for CS');

-- --------------------------------------------------------

--
-- Table structure for table `Grades`
--

CREATE TABLE IF NOT EXISTS `Grades` (
  `StudentID` int(10) unsigned NOT NULL,
  `Course` int(10) unsigned NOT NULL,
  `Grade` tinyint(4) DEFAULT NULL,
  `Term` varchar(20) COLLATE utf8_unicode_ci DEFAULT NULL
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci COMMENT='Grades in Python U.';

--
-- Dumping data for table `Grades`
--

INSERT INTO `Grades` (`StudentID`, `Course`, `Grade`, `Term`) VALUES
(1, 1, 7, '2016/1'),
(1, 2, 8, '2016/1'),
(2, 1, 6, '2015/1'),
(2, 2, 9, '2015/2'),
(3, 1, 5, '2016/1'),
(3, 2, 7, '2016/1'),
(4, 1, 9, '2016/1'),
(4, 2, 8, '2016/2');

-- --------------------------------------------------------

--
-- Table structure for table `Students`
--

CREATE TABLE IF NOT EXISTS `Students` (
  `ID` int(10) unsigned NOT NULL,
  `Name` varchar(150) COLLATE utf8_unicode_ci NOT NULL,
  `LastName` varchar(200) COLLATE utf8_unicode_ci NOT NULL,
  `DateJoined` date NOT NULL,
  `OutstandingBalance` tinyint(1) NOT NULL
) ENGINE=MyISAM AUTO_INCREMENT=6 DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci COMMENT='Main table.';

--
-- Dumping data for table `Students`
--

INSERT INTO `Students` (`ID`, `Name`, `LastName`, `DateJoined`, `OutstandingBalance`) VALUES
(1, 'Harry', 'Wilkinson', '2016-02-10', 0),
(2, 'Jonathan', 'Hunt', '2014-02-16', 0),
(3, 'Harry', 'Hughes', '2015-03-20', 0),
(4, 'Kayla', 'Allen', '2016-03-15', 1),
(5, 'Virginia', 'Gonzalez', '2017-04-02', 0);

--
-- Indexes for dumped tables
--

--
-- Indexes for table `Courses`
--
ALTER TABLE `Courses`
  ADD PRIMARY KEY (`CourseID`);

--
-- Indexes for table `Students`
--
ALTER TABLE `Students`
  ADD PRIMARY KEY (`ID`);

--
-- AUTO_INCREMENT for dumped tables
--

--
-- AUTO_INCREMENT for table `Courses`
--
ALTER TABLE `Courses`
  MODIFY `CourseID` int(10) unsigned NOT NULL AUTO_INCREMENT,AUTO_INCREMENT=3;
--
-- AUTO_INCREMENT for table `Students`
--
ALTER TABLE `Students`
  MODIFY `ID` int(10) unsigned NOT NULL AUTO_INCREMENT,AUTO_INCREMENT=6;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
