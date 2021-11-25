plugins {
    kotlin("jvm") version "1.6.0"
    `maven-publish` apply true
    `java-library` apply true
}

group = "com.github.agusotto96"
version = "1.0-SNAPSHOT"

repositories {
    mavenCentral()
}

dependencies {
    testImplementation("org.junit.jupiter:junit-jupiter:5.8.1")
}

tasks.test {
    useJUnitPlatform()
}